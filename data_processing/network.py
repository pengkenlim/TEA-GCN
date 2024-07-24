#setting sys.path for importing modules
import os
import sys
import scipy
import math
import numpy as np

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

from coexpression import ensemble

def calculate_edge_index(s, t, num_nodes):
    if s > t:
        s, t = t, s  # Ensure s < t
        return False, num_nodes * s - ((s + 1) * s) // 2 + t - s - 1
    else:
        return True, num_nodes * s - ((s + 1) * s) // 2 + t - s - 1

def init_flat_half_adj(num_edges):
    edge_values = np.empty(num_edges)
    edge_values[:] = math.nan
    return edge_values

def init_adj(num_genes):
    edge_values = np.empty( (num_genes, num_genes) )
    edge_values[:] = math.nan
    return edge_values


def extract_ranks(network_dir):
    genes = os.listdir(network_dir)
    gene_dict = {gene:i for i, gene in enumerate(genes)}
    num_genes = len(genes)
    num_edges = num_genes * (num_genes - 1) // 2
    upper_half = init_flat_half_adj(num_edges)
    lower_half = init_flat_half_adj(num_edges)
    upper_half_cor = init_flat_half_adj(num_edges)
    for idx, source in enumerate(genes):
        with open(os.path.join(network_dir, source), "r") as f:
            for line_no, line in enumerate(f):
                if line_no != 0 and line != "": #skip first and last line
                    target, cor, Rank = line.split("\n")[0].split("\t")
                    if target != source:
                        source_idx = gene_dict[source]
                        target_idx = gene_dict[target]
                        upper, edge_index= calculate_edge_index(source_idx, target_idx, num_genes)
                        if upper:
                            upper_half[edge_index] = Rank
                        else:
                            lower_half[edge_index] = Rank
                        upper_half_cor[edge_index] = cor
        if idx%1000 == 0:
            print(idx,"genes loaded")
    return genes , upper_half, lower_half, upper_half_cor, gene_dict

def write_all_edge_attributes(gene_dict, HRR_array, MR_array, cor_zscore_array, HRR_zscore_array, MR_zscore_array, old_network_dir, new_network_dir, genes):
    for idx, source in enumerate(genes):
        with open(os.path.join(old_network_dir, source), "r") as fin:
            with open(os.path.join(new_network_dir, source), "w") as fout:
                for line_no, line in enumerate(fin):
                    if line_no == 0:
                        fout.write("Target\tCo-exp_Str\tCo-exp_Str_ranking\tCo-exp_Str_HRR\tCo-exp_Str_MR\tzScore(Co-exp_Str)\tzScore(Co-exp_Str_HRR)\tzScore(Co-exp_Str_MR)\n")
                    elif line != "":
                        target, cor, Rank = line.split("\n")[0].split("\t")
                        cor = np.round(float(cor), 5)
                        if target == source:
                            MR , HRR , zHRR = 1.0 , 1.0, math.nan
                        else:
                            source_idx = gene_dict[source]
                            target_idx = gene_dict[target]
                            _, edge_index= calculate_edge_index(source_idx, target_idx, len(genes))
                            MR , HRR, zcor, zHRR, zMR = MR_array[edge_index] , HRR_array[edge_index], cor_zscore_array[edge_index], HRR_zscore_array[edge_index], MR_zscore_array[edge_index]
                            MR , HRR, zcor, zHRR , zMR = np.round(MR, 3) , np.round(HRR, 3), np.round(zcor, 7) , np.round(zHRR,7), np.round(zMR,7)
                        fout.write(f"{target}\t{cor}\t{Rank}\t{HRR}\t{MR}\t{zcor}\t{zHRR}\t{zMR}\n")
        if idx%1000 == 0:
            print(idx,"genes done")

def write_MR_attributes_only(gene_dict, MR_array,  MR_zscore_array, old_network_dir, new_network_dir, genes):
    for idx, source in enumerate(genes):
        with open(os.path.join(old_network_dir, source), "r") as fin:
            with open(os.path.join(new_network_dir, source), "w") as fout:
                for line_no, line in enumerate(fin):
                    if line_no == 0:
                        fout.write("Target\tCo-exp_Str\tCo-exp_Str_MR\tzScore(Co-exp_Str_MR)\n")
                    elif line != "":
                        target, cor, Rank = line.split("\n")[0].split("\t")
                        cor = np.round(float(cor), 5)
                        if target == source:
                            MR , HRR , zHRR = 1.0 , 1.0, math.nan
                        else:
                            source_idx = gene_dict[source]
                            target_idx = gene_dict[target]
                            _, edge_index= calculate_edge_index(source_idx, target_idx, len(genes))
                            MR , zMR = MR_array[edge_index] , MR_zscore_array[edge_index]
                            MR , zMR = np.round(MR, 3) , np.round(zMR,7)
                        fout.write(f"{target}\t{cor}\t{MR}\t{zMR}\n")
        if idx%1000 == 0:
            print(idx,"genes done")

