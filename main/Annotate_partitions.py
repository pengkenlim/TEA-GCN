#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse

# import argparse to parse thread information so that we can set thread environment variable before importing numpy and sklearn modules
if __name__ == "__main__":
    parser= argparse.ArgumentParser(description="Annotate_partitions.py\n\
                                        Lemmatize metadata of RNA-seq samples. Annotate dataset partitions with overrepresented lemmas")
        
    parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
    help = "Directory to output. Must be the same as for Generate_partitions.py, Run_TEA-GCN.py and Rank_transform.py.")


    parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
    help = "Delimiter for expression matrix. -de=\"t\" for tab-separated (.tsv). -de=\"c\" for comma separated (.csv). TSV by default." )

    parser.add_argument("-im", "--input_matrix_path", type=str, metavar="", required = True,
    help = "Path of expression matrix to input" )

    parser.add_argument("-k", "--k_clusters", type=int, metavar="", required = True,
    help = "Number of clusters to partition the expression data. For k=0, will use best k as determined in Generate_partitions.py step." )

    parser.add_argument("-m", "--model", type=str, metavar="", required = False, default ="en_core_web_sm" ,
    help = "spaCy model to use for lemmatization. '\'en_core_web_sm'\' will be used by default." )

    parser.add_argument("-p", "--permutations", type=int, metavar="", required = False, default =10_000 ,
    help = "Number of permutations to shuffle during calculation of overrepresentation p-values. Default = 10,000")

    args=parser.parse_args()
    output_dir = args.output_dir
    delimiter = args.delimiter
    input_matrix_path = args.input_matrix_path
    k_clusters = args.k_clusters
    model = args.model
    permutations = args.permutations

    import pickle
    import numpy as np
    from data_processing import nlp, expression_matrix

    sub_outdir =  os.path.join(output_dir, "Lemma_annotations")
    expmatfile = input_matrix_path

    delimiter = args.delimiter
    if delimiter == "t":
        delim = "\t"
    else:
        delim= ","

    if k_clusters ==0:
        k_clusters = read_write.load_pickle( os.path.join(output_dir,"Generate_partitions","selected_k.pkl"))
        print(f"k=0 specified. Will use best k of {k_clusters} as determined by Generate_partitions.py")
        k= int(k_clusters)
    else:
        k = k_clusters
    
    #load partition info
    with open(os.path.join(output_dir, "Generate_partitions", "k_cluster_assignment_dict.pkl") , "rb") as fin:
        k_cluster_assignment_dict = pickle.load(fin)
    cluster_assignment = k_cluster_assignment_dict[k]

    metadata_path = os.path.join(sub_outdir,"metadata.tsv")
    metadata_lemma_out_path  = os.path.join(sub_outdir,"metadata_lemmas.tsv")
    print("Proceeding to lemmatise metadata of RNA-seq samples...")
    nlp.lemmatise_metadata(metadata_path, metadata_lemma_out_path)
    print(f"Lemmatization complete. Lemmas of RNA-seq samples written to {metadata_lemma_out_path}")

    print("Proceeding to annotate dataset partitions with overrepresented lemmas...")
    
    #pre-counting
    sample_list = expression_matrix.load_only_samples(expmatfile, expmatsep = "\t")
    lemma_dict = nlp.load_metadata_lemma(metadata_lemma_out_path)
    sample_counts_dict, lemmas_index_dict = nlp.count_by_samples(lemma_dict)

    #determine observed counts
    partition_counts = nlp.generate_partition_counts(cluster_assignment, sample_list, sample_counts_dict)

    #first_shuffle
    shuffled_cluster_assignment = nlp.shuffle_cluster_assignment(cluster_assignment)
    P_value_matrix = np.ones((len(partition_counts), len(lemmas_index_dict)), dtype=int)
    
    #for every permutation...
    for i in range(permutations):
        #shuffle
        shuffled_cluster_assignment = nlp.shuffle_cluster_assignment(shuffled_cluster_assignment)
        #get shuffled counts
        shuffled_partition_counts = nlp.generate_partition_counts(shuffled_cluster_assignment, sample_list, sample_counts_dict)
        #update p_values
        P_value_matrix = nlp.compared_counts(partition_counts, shuffled_partition_counts, lemmas_index_dict, P_value_matrix)
        if i % 1000 ==0:
            print(f"{i} permutations completed.")
    P_value_matrix = P_value_matrix/ permutations

    
    P_value_matrix_path = os.path.join(sub_outdir, "P_value_matrix.pkl")
    lemmas_index_dict_path = os.path.join(sub_outdir, "lemmas_index_dict.pkl")
    with open( P_value_matrix_path, "wb") as fout:
         pickle.dump(P_value_matrix_path, fout)
    with open( lemmas_index_dict_path, "wb") as fout:
         pickle.dump(lemmas_index_dict_path, fout)

    partition_annotations_path = os.path.join(sub_outdir, "partition_annotation.tsv")
    all_OR_lemmas = nlp.write_annotation( P_value_matrix, lemmas_index_dict, partition_annotations_path)
    
    overrepresented_lemmas_path = os.path.join(sub_outdir, "overrepresented_lemmas.tsv")
    with open(overrepresented_lemmas_path, "w") as fout:
        for lemma in all_OR_lemmas:
            fout.write(lemma + "\n")
         
    print(f"Annotate_partitions.py complete.\n\
          Partition annotations written to {partition_annotations_path}\n\
        Overrepresented lemmas written to {overrepresented_lemmas_path}")


