#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse

if __name__ == "__main__":
        
    parser= argparse.ArgumentParser(description="Prepare_neighbourhoods.py. Prepare co-expression neighbourhoods of Genes-of-intrests from completed TEA-GCNs for gene function prediction using the GSEA algorithm.")
        
    parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
        help = "Directory to output. Must be the same as for Generate_partitions.py, Run_TEA-GCN.py and Rank_transform.py.")

    parser.add_argument("-g", "--gene_list", type=str, metavar="", required = True,
                             help = "list of comma-separated genes. spaces will be removed. E.g. \'GENE1,GENE2,GENE3\'")
    
    parser.add_argument("-n", "--net_name" , type=str, metavar = "" , default = "auto",
                        help = "name of the TEA-GCN to extract. Set to \"auto\" if you only generated one TEA-GCN. \"auto\" by default")
    
    args=parser.parse_args()
    output_dir=args.output_dir
    gene_list = args.gene_list
    net_name = args.net_name
    net_dir = os.path.join(output_dir, "Completed_GCN")

    from data_processing import read_write
    from scipy import stats
    import numpy as np

    if net_name == "auto":
        try:
            net_names = os.listdir(net_dir)
            if len(net_names) != 1:
                 print(f"Unable to accept net_name = \"auto\". More than one network or no network found. PLease check {net_dir} to make sure that there is only one directory.")
                 sys.exit()
            else:
                 net_path = os.path.join(net_dir, net_names[0])
        except:
             sys.exit(f"{net_dir} does not exist. Please check again.")
    else:
        #check net_name exists
        net_path = os.path.join(net_dir, net_name)
        if not os.path.exists(net_path):
             sys.exit(f"No network found in {net_path}. Please check again.")

    genes = gene_list.replace(" ","").split(",") #remove spaces and split
    temp_genes = []
    all_genes = os.listdir(net_path)
    for gene in genes:
        if gene not in all_genes:
             print(f"{gene} not found in {net_path}. Skipping {gene}...")
        else:
            temp_genes.append(gene)
    genes = temp_genes

    sub_outdir = os.path.join(output_dir, "Co-exp_Neighbourhoods")
    read_write.establish_dir(sub_outdir, isdir=True)

    import math

    for gene in genes:
        scores = []
        targets = []
        gene_path = os.path.join(net_path, gene)
        with open(gene_path, "r") as fin:
            for line_no, line in enumerate(fin):
                if line_no !=0 and line!= "" and line!= "\n":
                    score = float(line.split("\t")[-1])
                    target = (line.split("\t")[0])
                    if np.isnan(score):
                        score = 0
                    scores.append(score)
                    targets.append(target)

        ranks = stats.rankdata( np.array(scores) , method = "min", nan_policy = "omit")
        ranks_centered = ranks - int(np.nanmedian(ranks))

        new_gene_path = os.path.join(sub_outdir, gene)
        with open(new_gene_path, "w") as fout:
            fout.write("Target\tScore\tCentered_rank\n")
            for target, score, rank_centered in zip(targets, scores, ranks_centered):
                fout.write(f"{target}\t{score}\t{int(rank_centered)}\n")
        print(f"Co-expression neighbourhood of {gene} generated.")
                 


        