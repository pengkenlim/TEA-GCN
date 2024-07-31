#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

#parent_module = "/home/ken/Plant-GCN/src/" # remove
#sys.path.insert(0, parent_module) # remove

import numpy as np
import math
import scipy
from coexpression import ensemble , pearson , spearman, bicor
import warnings

def precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter="\t", workers=2):
    print("Precalculations for each correlation coefficient...")
    PCC_genes, PCC_gene_dict, PCC_nominators_dict, PCC_denominators_dict = pearson.precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delimiter, workers=workers)
    print("PCC precalculations complete")
    SCC_genes, SCC_gene_dict, SCC_nominators_dict, SCC_denominators_dict = spearman.precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delimiter, workers=workers)
    print("SCC precalculations complete")
    bicor_genes, bicor_gene_dict , bicor_norm_weights_dict = bicor.precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delimiter, workers=workers)
    print("bicor precalculations complete")
    return PCC_genes, PCC_gene_dict,SCC_genes, SCC_gene_dict, bicor_genes, bicor_gene_dict , PCC_nominators_dict, PCC_denominators_dict, SCC_nominators_dict, SCC_denominators_dict, bicor_norm_weights_dict

def calc_job(k, aggregation_method, PCC_genes,SCC_mapping_array,  bicor_mapping_array , SCC_mapping_dict, bicor_mapping_dict, network_path, gene_idx, gene, PCC_nominators_dict, PCC_denominators_dict,SCC_nominators_dict , SCC_denominators_dict,bicor_norm_weights_dict, threads ,full=False):
    warnings.filterwarnings(action='ignore', message='All-NaN slice encountered') 
    All_cor_means = []
    for cluster in range(k):
        PCC_cor_means = pearson.one_v_all(gene_idx, cluster, PCC_nominators_dict, PCC_denominators_dict)
        SCC_cor_means = spearman.one_v_all(SCC_mapping_dict[gene], cluster, SCC_nominators_dict, SCC_denominators_dict)[SCC_mapping_array]
        bicor_cor_means = bicor.one_v_all(bicor_mapping_dict[gene], cluster, bicor_norm_weights_dict, threads)[bicor_mapping_array]
        combined_cor_means = np.nanmax(np.array([PCC_cor_means , SCC_cor_means, bicor_cor_means]), axis =0 )
        All_cor_means.append(combined_cor_means)
    All_cor_means = np.array(All_cor_means)
    
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
        for target, cor, ES, rank in zip(PCC_genes, All_cor_means, ensemble_scores ,ensemble_ranks):
            if full:
                #f.write(f"{target}\t{cor}\t{ES}\t{rank}\n") #functionality not fleshed out yet
                f.write(f"{target}\t{ES}\t{rank}\n")
            else:
                f.write(f"{target}\t{ES}\t{rank}\n")
    return f"Calculated TEA co-expression strength ({aggregation_method}) for gene:{gene}, {gene_idx} out of {len(PCC_genes)}"

def array_mapper(ref_array ,  pertrubed_array_1, pertrubed_array_2):
    if type(ref_array) == list:
        ref_array = np.array(ref_array)
    if type(pertrubed_array_1) == list:
        pertrubed_array_1 = np.array(pertrubed_array_1)
    if type(pertrubed_array_1) == list:
        pertrubed_array_2 = np.array(pertrubed_array_2)
    unique_values_1 = {val: idx for idx, val in enumerate(pertrubed_array_1)}
    unique_values_2 = {val: idx for idx, val in enumerate(pertrubed_array_2)}
    mapping_array_1 =[]
    mapping_array_2 =[]
    for idx, val in enumerate(ref_array):
        mapping_array_1.append(unique_values_1[val])
        mapping_array_2.append(unique_values_2[val])
    mapping_array_1 = np.array(mapping_array_1)
    mapping_array_2 = np.array(mapping_array_2)
    return mapping_array_1,  mapping_array_2, unique_values_1, unique_values_2

def calc_untargeted(k, PCC_genes, SCC_genes, bicor_genes, PCC_nominators_dict, PCC_denominators_dict,SCC_nominators_dict , SCC_denominators_dict,bicor_norm_weights_dict, aggregation_method, network_path, workers= 2):
    threads = workers
    SCC_mapping_array , bicor_mapping_array ,SCC_mapping_dict  , bicor_mapping_dict= array_mapper(PCC_genes, SCC_genes, bicor_genes)
    for gene_idx, gene in enumerate(PCC_genes):
        if True:
            print(gene_idx)
            result=calc_job(k, aggregation_method, PCC_genes,SCC_mapping_array,  bicor_mapping_array , SCC_mapping_dict, bicor_mapping_dict, network_path, gene_idx, gene, PCC_nominators_dict, PCC_denominators_dict,SCC_nominators_dict , SCC_denominators_dict,bicor_norm_weights_dict, threads ,full=False)
            print(result)


def build_ensemble_GCN( Tid2Gid_dict, k_cluster_assignment_dict, expmat_path, k, network_path, aggregation_method, delim, workers):
    PCC_genes, PCC_gene_dict,SCC_genes, SCC_gene_dict, bicor_genes, bicor_gene_dict , PCC_nominators_dict, PCC_denominators_dict, SCC_nominators_dict, SCC_denominators_dict, bicor_norm_weights_dict = precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delim, workers=workers)
    print("Calculating and writing correlations...")
    calc_untargeted(k, PCC_genes, SCC_genes, bicor_genes, PCC_nominators_dict, PCC_denominators_dict,SCC_nominators_dict , SCC_denominators_dict,bicor_norm_weights_dict, aggregation_method, network_path, workers= workers)

def calc_targeted(k, path, edges ,PCC_genes, PCC_gene_dict,SCC_genes, SCC_gene_dict, bicor_genes, bicor_gene_dict ,PCC_nominators_dict, PCC_denominators_dict,SCC_nominators_dict , SCC_denominators_dict,bicor_norm_weights_dict,threads ,workers= 2):
    warnings.filterwarnings(action='ignore', message='All-NaN slice encountered') 
    #add header if creating new file
    if not os.path.exists(path):
        with open(path, "w") as f:
            f.write("Edge\tPartition co-exp.\n")
    PCC_source_array , PCC_target_array= [], []
    SCC_source_array , SCC_target_array= [], []
    bicor_source_array , bicor_target_array= [], []
    for edge in list(edges):
        source, target = edge.split("--")
        PCC_source_array.append(PCC_gene_dict[source])
        PCC_target_array.append(PCC_gene_dict[target])
        
        SCC_source_array.append(SCC_gene_dict[source])
        SCC_target_array.append(SCC_gene_dict[target])
        
        bicor_source_array.append(bicor_gene_dict[source])
        bicor_target_array.append(bicor_gene_dict[target])
    
    ALL_cor_means = []
    for cluster in range(k):
        PCC_cor_means = pearson.calc_job_k(PCC_source_array, PCC_target_array, PCC_nominators_dict, PCC_denominators_dict, cluster, threads)
        SCC_cor_means = spearman.calc_job_k(SCC_source_array, SCC_target_array, SCC_nominators_dict, SCC_denominators_dict, cluster, threads)
        bicor_cor_means = bicor.calc_job_k(bicor_source_array, bicor_target_array, bicor_norm_weights_dict, cluster, threads)
        combined_cor_means = np.nanmax( np.array( [ PCC_cor_means, SCC_cor_means , bicor_cor_means] ) , axis = 0)
        ALL_cor_means.append(combined_cor_means)

    ALL_cor_means = np.array(ALL_cor_means)
    
    #batch_ensemble_scores = []
    #for mode in ["Max", "Avg", "RAvg", "RWA", "RRWA"]:
        #batch_ensemble_scores.append(ensemble.aggregate(ALL_cor_means, mode, axis = 0))
    #batch_ensemble_scores = np.array(batch_ensemble_scores)
    with open(path, "a") as f:
        for idx, edge in enumerate(edges):
            cluster_cor = ALL_cor_means[:,idx]
            #ensemble_scores = batch_ensemble_scores[:,idx]
            cluster_cor = ",".join([str(i) for i in cluster_cor])
            #ensemble_scores = ",".join([str(i) for i in ensemble_scores])
            f.write(f"{edge}\t{cluster_cor}\n")

def get_partition_coexp(k, positive_met_edges_cor_path ,expmat_path, Tid2Gid_dict,  k_cluster_assignment_dict, delim, workers, positive_met_edges, threads):
    PCC_genes, PCC_gene_dict,SCC_genes, SCC_gene_dict, bicor_genes, bicor_gene_dict , PCC_nominators_dict, PCC_denominators_dict, SCC_nominators_dict, SCC_denominators_dict, bicor_norm_weights_dict = precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delim, workers=workers)
    print("Calculating and writing correlations of edges-of-interests...")

    if len(positive_met_edges) < 5000: 
        calc_targeted(k, positive_met_edges_cor_path, positive_met_edges ,
                    PCC_genes, PCC_gene_dict,SCC_genes, SCC_gene_dict, bicor_genes, bicor_gene_dict ,  PCC_nominators_dict, PCC_denominators_dict,SCC_nominators_dict , SCC_denominators_dict,bicor_norm_weights_dict,threads ,workers=workers)

    else:
        positive_met_edges_chunks = [positive_met_edges[x: x+ 5000] for x in range(0, len(positive_met_edges), 5000)]
        for chunk in positive_met_edges_chunks:
            calc_targeted(k, positive_met_edges_cor_path, chunk , 
                          PCC_genes, PCC_gene_dict,SCC_genes, SCC_gene_dict, bicor_genes, bicor_gene_dict ,  PCC_nominators_dict, PCC_denominators_dict,SCC_nominators_dict , SCC_denominators_dict,bicor_norm_weights_dict,threads ,workers=workers)

