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
        parser= argparse.ArgumentParser(description="Generate_partitions.py\n\
                                        Load expression matrix. Perform sample standardization to correct batch effect. Embed samples from gene-space into PC-space. Perform k-means clustering across a\
                                        range of K. Outputs principal component embeddings of samples (PCA data), mean silhouette coefficients from each clustering iteration, k-means cluster assignments of samples in pickle format.")

        parser.add_argument("-ks", "--k_start", type=int, metavar="", default=10, 
        help = "Starting k for K-means clustering. The default is 10")
        
        parser.add_argument("-ke", "--k_end", type=int, metavar="", default=50, 
        help = "Ending k for K-means clustering. The default is 50")

        parser.add_argument("-st", "--step", type=int, metavar="", default=1,
        help = "Step to increase between each k. The user is suggested to increase the step if the range is very big to decrease search space." )
        
        parser.add_argument("-t", "--threads", type=int, metavar="", default=4,
        help = "Number of threads to use for PCA and Kmeans clustering. Handled by Sklearn's parallelization. Try to specify lower than 64 for now." )

        parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
        help = "Working directory to output data." )

        parser.add_argument("-im", "--input_matrix_path", type=str, metavar="", required = True,
        help = "Path of expression matrix to input" )

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = " Delimiter for expression matrix. -de=\"t\" for tab-separated (.tsv). -de=\"c\" for comma separated (.csv). TSV by default." )
        
        parser.add_argument("-pc", "--pc_no", type= int, metavar="", default = 100,
        help = "Number of Principal components to keep. Default = 100. It cannot be more than the number of samples and the number of genes in the expression matrix." )

        parser.add_argument("-rs", "--randomstate", type= int, metavar="", default = 42,
        help = "Randomstate (seed) for random seeding during k-means clustering. By default, 42 will be used." )

        args=parser.parse_args()
        threads = args.threads
        
        #set threads and then import
        os.environ["MKL_NUM_THREADS"] = str(threads)
        os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
        os.environ["OMP_NUM_THREADS"] = str(threads)
        from sklearn_wrappers import kmeans , PCA
        from data_processing import read_write, expression_matrix
        
        k_start = args.k_start
        k_end = args.k_end
        step = args.step
        input_matrix_path = args.input_matrix_path
        pc_no = args.pc_no
        randomstate = args.randomstate
        k_list = [i for i in range(k_start,k_end+1,step)]

        #check validity of input then create outputdir if not already exists
        ouput_dir = args.output_dir
        #input_path = args.input_path
        sub_outdir =  os.path.join(ouput_dir, "Generate_partitions")

        read_write.establish_dir(sub_outdir , isdir= True)

        delimiter = args.delimiter
        if delimiter == "t":
                delim = "\t"
        else:
                delim= ","
        

        print("Loading expression matrix as pandas DataFrame...")
        expmat_path = input_matrix_path
        expmat_df = expression_matrix.load(expmat_path , expmatsep = delim, indexcolname="THIS_ARGUMENT_IS_DEPRECATED")

        print(f"Performing standardization of expression matrix and embedding samples in {pc_no} principal components...")
        pca_data , pc_variances = PCA.standardize_transform(expmat_df, n_pcs=pc_no)
        geneids = list(expmat_df.index)
        del expmat_df
        print(f"PCA complete. {pc_no} PCs representing {sum(pc_variances)}% of variance retained. Writing PCA data to output..." )
        read_write.to_pickle(pca_data , os.path.join(sub_outdir,"PCA_data.pkl"))

        print(f"Performing K-means clustering....")
        k_cluster_assignment_dict_path = os.path.join(sub_outdir,"k_cluster_assignment_dict.pkl")
        silhouette_coefficients_dict_path = os.path.join(sub_outdir, "silhouette_coefficients.pkl")
        k_cluster_assignment_dict , silhouette_coefficients, centroids_dict = kmeans.iterate_over_krange(pca_data,
                                                                                                         k_list , 
                                                                                                         k_cluster_assignment_dict_path , 
                                                                                                         silhouette_coefficients_dict_path, 
                                                                                                         mode = "kmeans",randomstate=randomstate)
        
        
        selected_k = k_list[silhouette_coefficients.index(max(silhouette_coefficients))]
        max_sc = max(silhouette_coefficients)

        print(f"k={selected_k} selected based on having the highest silhouette coefficient of {max_sc}\n")

        print(f"Writing k-means clustering data to output folder...")
        read_write.to_pickle(k_cluster_assignment_dict, os.path.join(sub_outdir,"k_cluster_assignment_dict.pkl"))
        read_write.to_pickle(silhouette_coefficients, os.path.join(sub_outdir, "silhouette_coefficients.pkl"))
        read_write.to_pickle(geneids, os.path.join(sub_outdir, "geneids.pkl"))
        
        sc_report_path = os.path.join(sub_outdir, "silhouette_coefficients_across_k.tsv")
        with open(sc_report_path, "w") as fout:
            fout.write("k\tsilhouette_coefficient\n")
            for k, sc in zip(k_list, silhouette_coefficients):
                   fout.write(f"{k}\t{sc}\n")
                
        read_write.to_pickle(selected_k, os.path.join(sub_outdir,"selected_k.pkl"))

        print("Generate_partitions.py complete.")


