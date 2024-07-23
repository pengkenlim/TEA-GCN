#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import numpy as np
import concurrent.futures as cf
import multiprocessing as mp
from coexpression import ensemble
import warnings
import scipy
import math
from einsumt import einsumt

def calc_job_k(source_array, target_array, norm_weights_dict, cluster, threads):
    warnings.filterwarnings(action='ignore', message='Mean of empty slice')
    #cor_values = np.einsum("ijk, ijk -> ij", norm_weights_dict[cluster][:,source_array,:], norm_weights_dict[cluster][:,target_array, :])
    cor_values = einsumt("ijk, ijk -> ij", np.take(norm_weights_dict[cluster], source_array , axis =1),
                         np.take(norm_weights_dict[cluster], target_array , axis =1), pool = threads)
    cor_means = np.nanmean(cor_values, axis=0)
    return cor_means


def get_norm_weights(cluster, all_values, assignment):
    np.seterr(all="ignore")
    cluster_values = all_values[assignment == cluster]
    subvalues_D = np.array([cluster_values]) # shape is (1, len(cluster_values)) # so that it can handle 2D input in the future
    values_minus_med_D = subvalues_D - np.nanmedian(subvalues_D, axis=1).reshape(-1,1)
    UI_D= (values_minus_med_D)/(9*np.nanmedian(np.abs(values_minus_med_D), axis=1).reshape(-1,1))
    Identity_D = np.where(abs(UI_D) < 1, 1,0)
    nom_D = (values_minus_med_D)*((1-UI_D**2)**2)* Identity_D
    norm_weights = nom_D/np.sqrt(np.sum(nom_D**2, axis =1)).reshape(-1,1)
    return cluster , norm_weights
                

def precalc_job(idx, line, delimiter, k,Tid2Gid_dict, assignment):
    norm_weights_gene_dict = {cluster:[] for cluster in range(k)}
    parts = line.rstrip().split(delimiter)
    all_values = [float(i) for i in parts[1:]]
    Transcript_id = parts[0]
    if Transcript_id in Tid2Gid_dict.keys():
        gene = Tid2Gid_dict[Transcript_id]
        #genes.append(gene)
        all_values=np.array(all_values)
        for cluster in range(k):
            cluster , norm_weights = get_norm_weights(cluster, all_values, assignment)
            norm_weights_gene_dict[cluster] = norm_weights
        if idx % 5000 ==0:
                    print("k=", k,":", idx , "genes prepared.")
        return gene, norm_weights_gene_dict
    return "ERROR", "ERROR"

def precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter="\t", workers=2):
    if k == 1:
        assignment = np.zeros(len(list(k_cluster_assignment_dict.values())[0]),dtype=int)
    else:
        assignment= k_cluster_assignment_dict[k]
    genes = []
    norm_weights_dict = {cluster:[] for cluster in range(k)}
    with open(expmat_path, 'r') as fin:
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            results = [executor.submit(precalc_job, idx, line, delimiter, k,Tid2Gid_dict, assignment) for idx, line in enumerate(fin) if idx != 0]
            for f in cf.as_completed(results):
                gene, norm_weights_gene_dict = f.result()
                if gene != "ERROR":
                    for cluster in range(k):
                        norm_weights_dict[cluster].append(norm_weights_gene_dict[cluster])
                    genes.append(gene)
    print("k=", k,":","Transposing normaized weights...")
    for cluster in range(k):
        norm_weights_dict[cluster]= np.transpose(np.array(norm_weights_dict[cluster]), (1,0,2))
        print(f"Cluster {cluster} transposed.")
    gene_dict={}
    for idx , gene in enumerate(genes, start=0):
        gene_dict[gene] = idx
    return genes, gene_dict , norm_weights_dict

def one_v_all(gene_idx, cluster, norm_weights_dict, threads):
    warnings.filterwarnings(action='ignore', message='Mean of empty slice')
    cor_values = einsumt("ijk, ijk -> ij", norm_weights_dict[cluster], np.broadcast_to(norm_weights_dict[cluster][:,[gene_idx], :], norm_weights_dict[cluster].shape),
                            pool = threads)
    cor_means = np.nanmean(cor_values, axis=0)
    return cor_means

def calc_job(k, aggregation_method, genes, network_path, gene_idx, gene, norm_weights_dict ,threads ,full=False):
    All_cor_means = []
    for cluster in range(k):
        cor_means  = one_v_all(gene_idx, cluster, norm_weights_dict, threads)
        All_cor_means.append(cor_means)
    ensemble_scores = ensemble.aggregate(All_cor_means, aggregation_method, axis = 0)
    ensemble_ranks =  scipy.stats.rankdata(ensemble_scores, method="min", nan_policy= "omit") # this will give maximum rank
    
    #flip ranks to reverse ranks
    if np.nanmax(ensemble_ranks) == 1: #happens when all ensemble scores are nan for what ever reason
        ensemble_ranks = [math.nan for i in range(len(ensemble_ranks))]
    else:
        ensemble_ranks = np.nanmax(ensemble_ranks) - ensemble_ranks +1
    
    All_cor_means = [",".join([str(i2) for i2 in i]) for i in  np.transpose(All_cor_means)] # list of comma separated raw correlations

    with open(os.path.join(network_path, gene), "w") as f:
        f.write(f"Target\t{aggregation_method}\tRank_of_target\n")
        for target, cor, ES, rank in zip(genes, All_cor_means, ensemble_scores ,ensemble_ranks):
            if full:
                #f.write(f"{target}\t{cor}\t{ES}\t{rank}\n") #functionality not fleshed out yet
                f.write(f"{target}\t{ES}\t{rank}\n")
            else:
                f.write(f"{target}\t{ES}\t{rank}\n")
    return f"Calculated bicor {aggregation_method} for sequence:{gene}, {gene_idx} out of {len(genes)}"


def calc_untargeted(k, genes , norm_weights_dict, aggregation_method, network_path, workers= 2):
    threads = workers
    for gene_idx, gene in enumerate(genes):
        result = calc_job( k, aggregation_method , genes, network_path, gene_idx, gene, norm_weights_dict, threads)
        print(result)


def build_ensemble_GCN( Tid2Gid_dict, k_cluster_assignment_dict, expmat_path, k, network_path, aggregation_method, delim, workers):
    genes, gene_dict , norm_weights_dict = precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delim, workers=workers)
    print("Calculating and writing correlations...")
    calc_untargeted(k, genes, norm_weights_dict, aggregation_method, network_path ,workers=workers)