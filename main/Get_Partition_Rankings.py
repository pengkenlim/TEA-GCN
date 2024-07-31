#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse


if __name__ == "__main__":
        
        parser= argparse.ArgumentParser(description="Get_Partition_Rankings.py.")
        
        parser.add_argument("-w", "--workers", type=int, metavar="", default=8,
        help = "Number of workers for parallelization. Affects precalc of all coefficients. But only affects the calc of bicor." )
        
        parser.add_argument("-t", "--threads", type=int, metavar="", default=8,
        help = "Number of threads for numpy linear algebra operations." )

        parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
        help = "Directory to output." )

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for expression matrix. -de=\"t\" for tab-separated (.tsv). -de=\"c\" for comma separated (.csv). TSV by default." )

        parser.add_argument("-k","--k_clusters", metavar="", type=int, default = 0,
        help = "Number of clusters to partition the expression data. For k=0, will use best k as determined in Generate_partitions.py step.")
        
        parser.add_argument("-im", "--input_matrix_path", type=str, metavar="", required = True,
        help = "Path of expression matrix to input" )

        parser.add_argument("-ep", "--edges_path", type=str, metavar="", required = True,
        help = "Path to list of edges to calculate.")

        parser.add_argument("-on", "--outname_prefix", type=str, metavar="", required = False, default = "Edges_of_interest",
        help = "prefix of to output calculations. ")

        

        args=parser.parse_args()
        
        workers=args.workers
        output_dir=args.output_dir
        delimiter=args.delimiter
        k_clusters = args.k_clusters
        threads =  args.threads
        input_matrix_path = args.input_matrix_path

        edges_path = args.edges_path
        outname_prefix = args.outname_prefix


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
        from scipy import stats

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
        
        edges = []
        with open(edges_path, "r") as fin:
            for line in fin:
                if line != "" and line != "\n":
                    line_contents = line.split("\n")[0].split(",")
                    edges.append("--".join(line_contents))

        #establish output sub directories
        sub_outdir = os.path.join(output_dir , "Partition_co-exp")

        network_dir_name = "_".join(["TEA",
                            f"{k}k"])
        
        network_path = os.path.join(sub_outdir, network_dir_name)
        edges_out_path = os.path.join(network_path, outname_prefix + "_coexp.tsv")
        

        
        read_write.establish_dir(network_path , isdir=True)
        print(f"Calculating Partition co-exp. of {len(edges)} edges with these specifications:\n\
              k= {k}\n\
              written to {network_path}")
        
        #loading required data
        geneids = read_write.load_pickle(os.path.join(output_dir , "Generate_partitions" ,"geneids.pkl"))
        Tid2Gid_dict = {geneid:geneid for geneid in geneids}
        k_cluster_assignment_dict =  read_write.load_pickle(os.path.join(output_dir, "Generate_partitions", "k_cluster_assignment_dict.pkl"))
        expmat_path = input_matrix_path

        #check of edges are geneids
        checked_edges = []
        for edge in edges:
            source, target = edge.split("--")
            try:
                catch = Tid2Gid_dict[source]
                catch = Tid2Gid_dict[target]
                checked_edges.append(edge)
            except:
                 print(f"WARNING. {edge} is not a valid edge!")

        if len(checked_edges) == 0:
            sys.exit()

        #start
        TEA.get_partition_coexp(k, edges_out_path ,expmat_path, Tid2Gid_dict,  k_cluster_assignment_dict, delim, workers, checked_edges, threads)

        #load
        edges_out_ranked_path =  os.path.join(network_path, outname_prefix + "_centered_ranks.tsv")
        with open(edges_out_path, "r") as fin:
             with open(edges_out_ranked_path, "w") as fout:
                  fout.write("Edge\tCentered_Ranks\n")
                  for line_no , line in enumerate(fin):
                        if line_no != 0 and line != "\n" and line != "" :
                            line_contents = line.split("\n")[0].split("\t")
                            Edge = line_contents[0]
                            scores = np.array([float(score) if score != "nan" else 0 for score in line_contents[1].split(",")])
                            ranks = stats.rankdata( np.array(scores) , method = "min", nan_policy = "omit")
                            ranks_centered = ranks - int(np.nanmedian(ranks))
                            rank_string = ",".join([str(rank) for rank in  list(ranks_centered)])
                            fout.write(f"{Edge}\t{rank_string}\n")

        print(f"Get_Partition_Rankings.py completed")
        print(f"Centered ranks of input edges generated at {edges_out_ranked_path}.")
