"""
Run_TEA-GCN_fast.py

Builds the completed TEA-GCN in one step, replacing
    Run_TEA-GCN.py -cc TEA -am RAvg   followed by   Rank_transform.py (default MR attributes)
with the same output in Completed_GCN/<NETWORK_NAME>/, but without writing and
re-reading the genes x genes intermediate network (RAW_GCN) as text and with
every pass parallelised over genes.

Requires Generate_partitions.py to have been run on the same --output_dir.

Same maths as the original scripts:
  per partition c: PCC, SCC and bicor of every gene pair, combined = nanmax(PCC, SCC, bicor)
  Co-exp_Str  = RAvg over partitions (nan -> 0, negatives -> 0, mean)
  rank_s(t)   = rank of t in s's list (1 = strongest, scipy rankdata "min", flipped)
  MR(s, t)    = sqrt(rank_s(t) * rank_t(s))
  zScore(MR)  = -(MR - mean) / std over all gene pairs; self row MR 1.0 and the max zScore
Co-exp_Str is computed with the same numpy calls as TEA.calc_job, so it is
bit-identical to Run_TEA-GCN.py. This matters because small partitions produce
many mathematically tied scores whose order (and so their ranks and MR) is
decided by the last floating-point bit. The only expected difference is in the
last digits of the unrounded self-row zScore: the original's mean/std summation
order follows os.listdir() order of RAW_GCN.

Memory: about 4 * genes^2 bytes (int32 rank matrix) plus small per-process buffers.
Uses fork-based multiprocessing (Linux / macOS).
"""

import os
import sys

if __name__ == "__main__":
    # one BLAS thread per process: parallelism comes from the worker processes, and the
    # single-thread kernels keep the recomputed scores of pass 3 bit-identical to pass 1
    for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[_v] = "1"

import argparse
import math
import mmap
import multiprocessing as mp
import pickle
import shutil
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import scipy.stats

HEADER = "Target\tCo-exp_Str\tCo-exp_Str_MR\tzScore(Co-exp_Str_MR)\n"

# shared with the forked worker processes (read-only, except the rank matrix in pass 1)
G = {}


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def load_matrix(path, delim):
    """Same parsing as the precalc functions: line.rstrip().split(delim), float() per value."""
    genes, rows = [], []
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split(delim)
            genes.append(parts[0])
            rows.append([float(v) for v in parts[1:]])
    return genes, np.array(rows, dtype=np.float64)


def precalc(X, assignment, k):
    """Row-wise versions of pearson.precalc, spearman.precalc and bicor.get_norm_weights."""
    out = []
    with np.errstate(all="ignore"):
        for c in range(k):
            V = np.ascontiguousarray(X[:, assignment == c])
            n_c = V.shape[1]
            p_nom = V - (np.sum(V, axis=1) / n_c)[:, None]
            p_den = np.sqrt(np.sum(p_nom ** 2, axis=1))
            Rk = scipy.stats.rankdata(V, method="average", nan_policy="omit", axis=1)
            s_nom = Rk - (np.sum(Rk, axis=1) / n_c)[:, None]
            s_den = np.sqrt(np.sum(s_nom ** 2, axis=1))
            med = np.nanmedian(V, axis=1).reshape(-1, 1)
            vmm = V - med
            UI = vmm / (9 * np.nanmedian(np.abs(vmm), axis=1).reshape(-1, 1))
            ident = np.where(abs(UI) < 1, 1, 0)
            nom = vmm * ((1 - UI ** 2) ** 2) * ident
            W = nom / np.sqrt(np.sum(nom ** 2, axis=1)).reshape(-1, 1)
            # bicor weights in the (1, genes, samples) layout built by bicor.precalc
            out.append(dict(p_nom=p_nom, p_den=p_den, s_nom=s_nom, s_den=s_den,
                            W3=np.ascontiguousarray(W)[None, :, :]))
    return out


# ---------------------------------------------------------------- co-expression strength

def coexp_row(g):
    """Co-exp_Str of gene g against all genes; mirrors TEA.calc_job call for call
    (np.dot per gene for PCC/SCC, einsum for bicor, nanmax, ensemble.aggregate RAvg)."""
    rows = []
    for d in G["pre"]:
        pcc = np.dot(d["p_nom"], d["p_nom"][g]) / np.dot(d["p_den"], d["p_den"][g])
        scc = np.dot(d["s_nom"], d["s_nom"][g]) / np.dot(d["s_den"], d["s_den"][g])
        W3 = d["W3"]
        bic = np.nanmean(np.einsum("ijk, ijk -> ij", W3, np.broadcast_to(W3[:, [g], :], W3.shape)), axis=0)
        rows.append(np.nanmax(np.array([pcc, scc, bic]), axis=0))
    cor = np.array(rows)
    cor[np.isnan(cor)] = 0
    cor[cor < 0] = 0
    return np.nanmean(cor, axis=0)


def coexp_block(i0, i1):
    with np.errstate(all="ignore"):
        return np.array([coexp_row(g) for g in range(i0, i1)])


# ---------------------------------------------------------------- pass 1: ranks

def rank_block(i0, i1):
    es = coexp_block(i0, i1)
    r = scipy.stats.rankdata(es, method="min", axis=1)
    mx = r.max(axis=1, keepdims=True)
    G["R"][i0:i1] = mx - r + 1
    bad = (mx[:, 0] == 1)  # all scores tied: Run_TEA-GCN.py writes nan ranks
    if bad.any():
        G["R"][i0:i1][bad] = 0
    return i1 - i0


def rank_blocks(blocks):
    return sum(rank_block(i0, i1) for i0, i1 in blocks)


# ---------------------------------------------------------------- pass 2: MR mean / std

def mr_block(i0, i1):
    R = G["R"]
    with np.errstate(invalid="ignore"):
        mr = np.sqrt(R[i0:i1].astype(np.float64) * R[:, i0:i1].T.astype(np.float64))
    mr[mr == 0] = np.nan  # rank 0 marks nan ranks
    return mr


def upper_pairs(i0, i1):
    mr = mr_block(i0, i1)
    u = mr[np.arange(mr.shape[1])[None, :] > np.arange(i0, i1)[:, None]]
    return u[~np.isnan(u)]


def stats_sum(i0, i1):
    u = upper_pairs(i0, i1)
    return np.sum(u), u.size, u.min() if u.size else np.inf


def stats_sq(i0, i1):
    u = upper_pairs(i0, i1)
    return np.sum((u - G["mean"]) ** 2)


# ---------------------------------------------------------------- pass 3: write Completed_GCN

def write_blocks(blocks):
    genes, out_dir, min_z = G["genes"], G["out_dir"], G["min_z"]
    mean, std, zmax = G["mean"], G["std"], G["zmax"]
    zmax_s = repr(float(zmax))
    gene_arr = np.array(genes, dtype=object)
    n_rows = 0
    for i0, i1 in blocks:
        es = np.round(coexp_block(i0, i1), 5)
        mr = mr_block(i0, i1)
        with np.errstate(invalid="ignore"):
            z = -1 * ((mr - mean) / std)
        mr = np.round(mr, 3)
        z = np.round(z, 7)
        for r, s in enumerate(range(i0, i1)):
            if min_z is None:
                idx = np.arange(len(genes))
            else:
                keep = z[r] >= min_z
                keep[s] = zmax >= min_z
                idx = np.flatnonzero(keep)
            lines = [HEADER]
            for t, cv, mv, zv in zip(gene_arr[idx].tolist(), es[r, idx].tolist(),
                                     mr[r, idx].tolist(), z[r, idx].tolist()):
                if t == genes[s]:
                    lines.append(f"{t}\t{cv!r}\t1.0\t{zmax_s}\n")
                else:
                    lines.append(f"{t}\t{cv!r}\t{mv!r}\t{zv!r}\n")
            tmp = os.path.join(out_dir, genes[s] + ".writing")
            with open(tmp, "w") as f:
                f.write("".join(lines))
            os.replace(tmp, os.path.join(out_dir, genes[s]))
            n_rows += len(lines) - 1
    return n_rows


# ---------------------------------------------------------------- main

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    parser = argparse.ArgumentParser(description="Run_TEA-GCN_fast.py. Builds the completed TEA-GCN (Co-exp_Str, MR and zScore(MR)) "
                                     "in one step, equivalent to Run_TEA-GCN.py -cc TEA -am RAvg followed by Rank_transform.py. "
                                     "Requires prior data partitioning using Generate_partitions.py in order to run.")
    parser.add_argument("-o", "--output_dir", type=str, metavar="", required=True,
                        help="Directory to output. Must be the same as for Generate_partitions.py")
    parser.add_argument("-im", "--input_matrix_path", type=str, metavar="", required=True,
                        help="Path of expression matrix to input")
    parser.add_argument("-de", "--delimiter", type=str, metavar="", default="t", choices=["t", "c"],
                        help="Delimiter for expression matrix. -de=\"t\" for tab-separated (.tsv). -de=\"c\" for comma separated (.csv). TSV by default.")
    parser.add_argument("-k", "--k_clusters", type=int, metavar="", default=0,
                        help="Number of clusters to partition the expression data. For k=0, will use best k as determined in Generate_partitions.py step.")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=os.cpu_count(),
                        help="Number of worker processes. All CPUs by default.")
    parser.add_argument("-mz", "--min_zscore", type=float, metavar="", default=None,
                        help="Only write edges with zScore(Co-exp_Str_MR) >= this value (e.g. 1). All edges by default, as Rank_transform.py.")
    parser.add_argument("-b", "--block", type=int, metavar="", default=0,
                        help="Genes per work unit (0 = automatic).")
    args = parser.parse_args()

    t0 = time.time()
    gp = os.path.join(args.output_dir, "Generate_partitions")
    with open(os.path.join(gp, "k_cluster_assignment_dict.pkl"), "rb") as f:
        k_assign = pickle.load(f)
    k = args.k_clusters
    if k == 0:
        with open(os.path.join(gp, "selected_k.pkl"), "rb") as f:
            k = int(pickle.load(f))
        print(f"k=0 specified. Will use best k of {k} as determined by Generate_partitions.py")
    assignment = np.zeros(len(next(iter(k_assign.values()))), dtype=int) if k == 1 else k_assign[k]

    genes, X = load_matrix(args.input_matrix_path, "\t" if args.delimiter == "t" else ",")
    n = len(genes)
    if len(set(genes)) != n:
        sys.exit("Duplicate gene IDs in the expression matrix.")
    log(f"Loaded {n} genes x {X.shape[1]} samples, k={k}, {args.threads} processes")

    G["k"] = k
    G["pre"] = precalc(X, assignment, k)
    del X
    G["genes"] = genes
    G["min_z"] = args.min_zscore
    # shared anonymous mapping: forked pass-1 workers write rank rows the parent can read
    G["R"] = np.frombuffer(mmap.mmap(-1, max(n * n * 4, 1)), dtype=np.int32, count=n * n).reshape(n, n)
    log(f"Precalculation complete. Rank matrix: {n * n * 4 / 1e9:.1f} GB")

    b = args.block or int(max(8, min(256, 256e6 / (n * 8 * 10))))
    blocks = [(i, min(i + b, n)) for i in range(0, n, b)]
    chunks = [blocks[i::args.threads * 4] for i in range(args.threads * 4)]
    ctx = mp.get_context("fork")

    # pass 1: co-expression strength and ranks of every gene's neighbourhood
    with ctx.Pool(args.threads) as pool:
        done, step = 0, max(n // 20, 1)
        for m in pool.imap_unordered(rank_blocks, chunks):
            prev, done = done, done + m
            if done // step > prev // step:
                log(f"Co-expression strengths and ranks: {done}/{n} genes")

    # pass 2: mean / std of MR over all gene pairs
    with ThreadPoolExecutor(args.threads) as ex:
        parts = list(ex.map(lambda ij: stats_sum(*ij), blocks))
    count = sum(p[1] for p in parts)
    G["mean"] = math.fsum(p[0] for p in parts) / count
    with ThreadPoolExecutor(args.threads) as ex:
        G["std"] = math.sqrt(math.fsum(ex.map(lambda ij: stats_sq(*ij), blocks)) / count)
    G["zmax"] = -1 * ((min(p[2] for p in parts) - G["mean"]) / G["std"])
    log(f"MR of {count} edges: mean={G['mean']:.6f} std={G['std']:.6f}")

    # pass 3: write one file per source gene
    net = f"TEA_{k}k_RAvg"
    out_dir = os.path.join(args.output_dir, "Completed_GCN", net)
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    G["out_dir"] = out_dir
    with ctx.Pool(args.threads) as pool:
        rows = sum(pool.map(write_blocks, chunks))
    minutes = (time.time() - t0) / 60
    with open(os.path.join(args.output_dir, "Completed_GCN", f"{net}.summary.tsv"), "w") as f:
        f.write(f"genes\t{n}\nk\t{k}\nedges_written\t{rows}\nmin_zscore\t{args.min_zscore}\n"
                f"mr_mean\t{G['mean']!r}\nmr_std\t{G['std']!r}\nmr_pairs\t{count}\nminutes\t{minutes:.2f}\n")
    print(f"Run_TEA-GCN_fast.py finished running ({minutes:.1f} min). Completed TEA-GCN generated at {out_dir}")
