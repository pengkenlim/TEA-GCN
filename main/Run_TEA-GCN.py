#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse


if __name__ == "__main__":
        
        parser= argparse.ArgumentParser(description="Run_TEA-GCN.py.\n\
                                        Construct TEA-GCN from expression data (in the form of an expression matrix). Requires prior data partitioning using Generate_partitions.py in order to run.")
        
        parser.add_argument("-w", "--workers", type=int, metavar="", default=8,
        help = "Number of workers for parallelization. Affects precalc of all coefficients. But only affects calc of bicor." )
        
        parser.add_argument("-t", "--threads", type=int, metavar="", default=8,
        help = "Number of threads for numpy linear algebra operations. Affects PCC and SCC calc." )

        parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
        help = "Directory to output. Must be the same as for Generate_partitions.py" )

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for expression matrix. -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default." )

        parser.add_argument("-cc","--correlation_coefficient", type= str, metavar="", default = "TEA", choices=["PCC","SCC","bicor", "TEA"] ,
        help = "select \'TEA\' to enable Coefficient Aggregation. Select \"PCC\",\"SCC\",\"bicor\" for constructing Partition Aggregation-only GCN" )

        parser.add_argument("-am","--aggregation_method", type=str,  metavar="" , default = "RAvg", choices = ["Max","Avg","RAvg","RWA","RRWA"] ,
        help = "Default is RAvg.")

        parser.add_argument("-k","--k_clusters", metavar="", type=int, default = 0,
        help = "Number of clusters to partition the expression data. For k=0, will use best k as determined in Generate_partitions.py step")
        
        parser.add_argument("-im", "--input_matrix_path", type=str, metavar="", required = True,
        help = "Path of expression matrix to input" )

        args=parser.parse_args()
        
        workers=args.workers
        output_dir=args.output_dir
        delimiter=args.delimiter
        correlation_coefficient = args.correlation_coefficient
        aggregation_method = args.aggregation_method
        k_clusters = args.k_clusters
        threads =  args.threads
        input_matrix_path = args.input_matrix_path

        #set threads and then import
        os.environ["MKL_NUM_THREADS"] = str(threads)
        os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
        os.environ["OMP_NUM_THREADS"] = str(threads)

        from coexpression import bicor , pearson , spearman , TEA
        #from data_processing import read_write, network
        from data_processing import read_write
        #from analyses import network_performance
        import numpy as np
        import pandas as pd

        if delimiter == "t":
            delim = "\t"
        else:
            delim = ","
        
        
        if k_clusters ==0:
            k_clusters = read_write.load_pickle( os.path.join(output_dir,"Generate_partitions","selected_k.pkl"))
            print(f"k=0 specified. Will use best k of {k_clusters} as determined by Generate_partitions.py")
            k= int(k_clusters)
        else:
            k = k_clusters
        
                #establish output sub directories
        sub_outdir = os.path.join(output_dir , "RAW_GCN")
        cc_sub_outdir = os.path.join(sub_outdir, correlation_coefficient)

        network_dir_name = "_".join([correlation_coefficient,
                            f"{k}k",
                            aggregation_method])
        
        network_path = os.path.join(cc_sub_outdir, network_dir_name)

        
        read_write.establish_dir(network_path , isdir=True)
        print(f"Building ensemble network with these specifications:\n\n\
              Correlation coefficient= {correlation_coefficient}\n\
              k= {k}\n\
              Aggregation method= {aggregation_method}\n\n\
              written to {network_path}")
        
        #loading required data
        geneids = read_write.load_pickle(os.path.join(output_dir , "Generate_partitions" ,"geneids.pkl"))
        Tid2Gid_dict = {geneid:geneid for geneid in geneids}
        k_cluster_assignment_dict =  read_write.load_pickle(os.path.join(output_dir, "Generate_partitions", "k_cluster_assignment_dict.pkl"))
        expmat_path = input_matrix_path

        #start
        if correlation_coefficient == "bicor":
            build_ensemble_GCN = bicor.build_ensemble_GCN
        elif correlation_coefficient== "PCC":
            build_ensemble_GCN = pearson.build_ensemble_GCN
        elif correlation_coefficient == "SCC":
            build_ensemble_GCN = spearman.build_ensemble_GCN
        elif correlation_coefficient == "TEA":
            build_ensemble_GCN = TEA.build_ensemble_GCN

        build_ensemble_GCN( Tid2Gid_dict, k_cluster_assignment_dict, expmat_path, k, network_path, aggregation_method, delim, workers)
        