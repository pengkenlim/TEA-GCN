#setting sys.path for importing modules
import os
import sys
import scipy

import spacy
from collections import Counter
import pandas as pd
import math
import re
from spacy.matcher import PhraseMatcher
from nltk.collocations import BigramCollocationFinder, BigramAssocMeasures
import inflect
import string
from sklearn.feature_extraction.text import TfidfVectorizer
import time
import numpy as np
import statsmodels.stats.multitest
import pickle
import csv

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

from data_processing import expression_matrix

def shuffle_cluster_assignment(cluster_assignment, seed = 42):
    np.random.seed(seed)
    shuffled_arr = np.random.permutation(cluster_assignment)
    return shuffled_arr



def clean_text(text, digit_replacement_table, punc_replacement_table):
  return text.strip().translate(digit_replacement_table).translate(punc_replacement_table)

def lemmatise_metadata(metadata_path, metadata_lemma_out_path, model = "en_core_web_sm"):
    nlp = spacy.load(model)
    digit_replacement_table = str.maketrans('', '', string.digits)
    punc_replacement_table = str.maketrans({char: ' ' for char in string.punctuation})
    with open(metadata_path, "r") as fin:
         with open(metadata_lemma_out_path, "w") as fout:
            fout.write("Accession\tLemmas\n")
            for line_no, line in enumerate(fin):
                if line_no != 0 and line != "" and line != "\n":
                    lemmas = []
                    line_contents = line.split("\n")[0].split("\t")
                    sample_name = line_contents[0]
                    line_contents = line_contents[1:]
                    text = " ".join(line_contents).lower()
                    c_text = clean_text(text, digit_replacement_table, punc_replacement_table)
                    doc = nlp(c_text)
                    for token in doc:
                        if not token.is_stop and not token.is_space:
                            lemmas.append(token.lemma_ )
                    lemmas = ",".join(lemmas)
                    fout.write(f"{sample_name}\t{lemmas}\n")
                if line_no % 1000 == 0:
                    print(f"Metadata of {line_no} RNA-seq samples lemmatised.")  

def load_metadata_lemma(metadata_lemma_out_path):
    lemma_dict = {}
    with open(metadata_lemma_out_path, "r") as fin:
        for line_no , line in enumerate(fin):
            if line_no != 0:
                line_contents = line.split("\n")[0].split("\t")
                lemma_dict[line_contents[0]] = line_contents[1].split(",")
    return lemma_dict

def count_by_samples(lemma_dict):
    """6 mins / 5000 samples"""
    all_lemmas = []
    sample_counts_dict ={} 
    for sample, lemmas in lemma_dict.items():
        sample_counts_dict[sample] = Counter(lemmas)
        all_lemmas.extend(lemmas)

    all_lemmas = list(set(all_lemmas))
    lemmas_index_dict = { lemma: idx for idx, lemma in enumerate(all_lemmas) }
    return sample_counts_dict, lemmas_index_dict

def generate_partition_counts(cluster_assignment, sample_list, sample_counts_dict):
    partition_counts = {}
    for partition , sample in zip(list(cluster_assignment), sample_list):
        try:
            sample_counts = sample_counts_dict[sample]
        except:
            sample_counts = Counter()
        try:
            partition_counts[partition].update(sample_counts.copy())
        except:
            partition_counts[partition] = sample_counts.copy()
    return partition_counts


    
def generate_partition_lemmas(cluster_assignment, sample_list, lemma_dict, for_shuffle=False):
    all_lemmas = []
    partition_lemmas = {}
    for cluster , sample in zip(list(cluster_assignment), sample_list):
        try:
            sample_lemmas = lemma_dict[sample]
        except:
            sample_lemmas = []
        try:
            partition_lemmas[cluster].extend(sample_lemmas)
        except:
            partition_lemmas[cluster]= sample_lemmas
        all_lemmas.extend(sample_lemmas)
    if not for_shuffle:
        all_lemmas = list(set(all_lemmas))
        lemmas_index_dict = { lemma: idx for idx, lemma in enumerate(all_lemmas) }
        return partition_lemmas , lemmas_index_dict
    else:
        return partition_lemmas


def compared_counts(partition_lemma_counts, shuffled_partition_lemma_counts, lemmas_index_dict, P_value_matrix):
    P_value_matrix += 1
    for partition, lemma_counts in partition_lemma_counts.items():
        shuffled_lemma_counts = shuffled_partition_lemma_counts[partition]
        count_difference = lemma_counts - shuffled_lemma_counts
        for lemma in count_difference.keys():
            lemma_idx = lemmas_index_dict[lemma]
            P_value_matrix[partition, lemma_idx] -= 1
    return P_value_matrix
        

def write_annotation(P_value_matrix, lemmas_index_dict, outpath):
    all_OR_lemmas = []
    index_lemma_dict = {value: key for key, value in lemmas_index_dict.items()}
    lemma_array=[]
    for i in range(len(index_lemma_dict)):
        lemma_array.append(index_lemma_dict[i])
    lemma_array = np.array(lemma_array)
    with open(outpath, "w") as fout:
        fout.write(f"Dataset Partition\tOverrepresented lemmas\t")
        for partition in range(P_value_matrix.shape[0]):
            p_values = P_value_matrix[partition,:]
            indices = np.where(p_values <=0.05 )
            extracted_lemma_array = list(lemma_array[indices])
            extracted_p_values = list(p_values[indices])
            col_string = [f"({lemma} : {np.round(p_val , decimals = 4)})" for lemma, p_val in zip(extracted_lemma_array, extracted_p_values)]
            col_string = ",".join(col_string)
            fout.write(f"{partition}\t{col_string}\t")
            all_OR_lemmas.extend(extracted_lemma_array)
    return list(set(all_OR_lemmas))

        

#testing code
if __name__ == "__main__":
    k = 93
    sep = "\t"
    path = "/mnt/md0/ken/correlation_networks/taxid3702_5k/Generate_partitions/k_cluster_assignment_dict.pkl"
    expmatfile = "/mnt/md0/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5k/QC_expression_data/expression_matrix.tsv"
    with open(path, "rb") as fin:
        k_cluster_assignment_dict = pickle.load(fin)
    cluster_assignment = k_cluster_assignment_dict[k]

    metadata_path = "/mnt/md0/ken/correlation_networks/taxid3702_5k/Lemma_annotations/metadata.tsv"
    metadata_lemma_out_path  = "/mnt/md0/ken/correlation_networks/taxid3702_5k/Lemma_annotations/metadata_lemmas.tsv"
    
    lemmatise_metadata(metadata_path, metadata_lemma_out_path)
    lemma_dict = load_metadata_lemma(metadata_lemma_out_path)
    sample_list = expression_matrix.load_only_samples(expmatfile, expmatsep = "\t")

    sample_counts_dict, lemmas_index_dict = count_by_samples(lemma_dict)


    partition_counts = generate_partition_counts(cluster_assignment, sample_list, sample_counts_dict)

    #first_shuffle
    shuffled_cluster_assignment = shuffle_cluster_assignment(cluster_assignment)
    P_value_matrix = np.ones((len(partition_counts), len(lemmas_index_dict)), dtype=int)
    for i in range(10_000):
        shuffled_cluster_assignment = shuffle_cluster_assignment(shuffled_cluster_assignment)
        shuffled_partition_counts = generate_partition_counts(shuffled_cluster_assignment, sample_list, sample_counts_dict)

        P_value_matrix = compared_counts(partition_counts, shuffled_partition_counts, lemmas_index_dict, P_value_matrix)
        if i % 100 ==0:
            print(f"{i} permutations completed.")
    P_value_matrix = P_value_matrix / 10_000

    P_value_matrix[:, lemmas_index_dict["anther"] ]
    P_value_matrix[:, lemmas_index_dict["flower"] ]
    np.min(P_value_matrix[1, : ])

    
    #P_bh_matrix = []
    #for i in range(len(partition_counts)):
    #    rejected, corrected_p_values,_,_= statsmodels.stats.multitest.multipletests(P_value_matrix[i,:], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    #    P_bh_matrix.append(corrected_p_values)
    #P_bh_matrix = np.array(P_bh_matrix)

    

    



                
                                
                                 
                                

                             


