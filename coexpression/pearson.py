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
from einsumt import einsumt
from coexpression import ensemble
import warnings
import scipy
import math


def precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter="\t", workers=2):
    if k == 1:
        assignment = np.zeros(len(list( k_cluster_assignment_dict.values() ) [0] ) , dtype=int)
    else:
        assignment= k_cluster_assignment_dict[k]
    genes = []
    nominators_dict, denominators_dict, = {cluster:[] for cluster in range(k)}, {cluster:[] for cluster in range(k)}
    with open(expmat_path, 'r') as fin:
        for idx, line in enumerate(fin):
            if idx != 0:
                parts = line.rstrip().split(delimiter)
                all_values = np.array([float(i) for i in parts[1:]])
                Transcript_id = parts[0]
                if Transcript_id in Tid2Gid_dict.keys():
                    gene = Tid2Gid_dict[Transcript_id]
                    genes.append(gene)               
                    for cluster in range(k): # precalc job
                        cluster_values = all_values[assignment == cluster]
                        nomi = cluster_values - np.array([(np.sum(cluster_values)/len(cluster_values))])
                        denomi = np.sqrt(np.sum(nomi**2))
                        nominators_dict[cluster].append(nomi)
                        denominators_dict[cluster].append(denomi)
            if idx % 5000 == 0:
                print("k=", k,":", idx , "genes prepared.")
    gene_dict={}
    for idx , gene in enumerate(genes, start=0):
        gene_dict[gene] = idx
    for cluster in range(k):
        nominators_dict[cluster] = np.array(nominators_dict[cluster])
        denominators_dict[cluster] = np.array(denominators_dict[cluster])
    
    return genes, gene_dict, nominators_dict, denominators_dict



def one_v_all(gene_idx, cluster, nominators_dict, denominators_dict):
    warnings.filterwarnings(action='ignore', message='invalid value encountered in divide')
    nominator = np.dot(nominators_dict[cluster], nominators_dict[cluster][gene_idx])
    denominator = np.dot(denominators_dict[cluster], denominators_dict[cluster][gene_idx])
    cor_means = nominator/denominator
    return cor_means

def calc_job(k, aggregation_method, genes, network_path, gene_idx, gene,   nominators_dict, denominators_dict ,full=False):
    All_cor_means = []
    for cluster in range(k):
        cor_means  = one_v_all(gene_idx, cluster, nominators_dict, denominators_dict)
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
    return f"Calculated PCC {aggregation_method} for sequence:{gene}, {gene_idx} out of {len(genes)}"

def calc_untargeted(k, genes, nominators_dict, denominators_dict, aggregation_method, network_path):
    for  gene_idx, gene in enumerate(genes):
        result=calc_job(k, aggregation_method , genes, network_path, gene_idx, gene, nominators_dict, denominators_dict)
        print(result)

def build_ensemble_GCN( Tid2Gid_dict, k_cluster_assignment_dict, expmat_path, k, network_path, aggregation_method, delim, workers):
    genes, gene_dict, nominators_dict, denominators_dict = precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delim, workers=workers)
    print("Calculating and writing correlations...")
    calc_untargeted(k, genes, nominators_dict, denominators_dict, aggregation_method, network_path)

def calc_job_k(source_array, target_array, shared_nominators_dict, shared_denominators_dict, cluster, threads):
    warnings.filterwarnings(action='ignore', message='invalid value encountered in divide')
    #cor_values = np.sum(shared_nominators_dict[cluster][source_array] * shared_nominators_dict[cluster][target_array], axis=1)/(shared_denominators_dict[cluster][source_array]* shared_denominators_dict[cluster][target_array])
    numerator = einsumt('ij,ij->i', np.take(shared_nominators_dict[cluster], source_array , axis = 0) , np.take(shared_nominators_dict[cluster], target_array , axis = 0), pool =threads)
    denominator = np.einsum('i,i->i',  np.take(shared_denominators_dict[cluster], source_array), np.take(shared_denominators_dict[cluster], target_array))
    cor_values = numerator / denominator
    return cor_values