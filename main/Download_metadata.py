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
    parser= argparse.ArgumentParser(description="Download_metadata.py\n\
                                        Download metadata of SRR accessions from European Nucleotide Archive (ENA).")
        
    parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
    help = "Directory to output. Must be the same as for Generate_partitions.py, Run_TEA-GCN.py and Rank_transform.py.")


    parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
    help = "Delimiter for expression matrix. -de=\"t\" for tab-separated (.tsv). -de=\"c\" for comma separated (.csv). TSV by default." )

    parser.add_argument("-im", "--input_matrix_path", type=str, metavar="", required = True,
    help = "Path of expression matrix to input" )

    
    args=parser.parse_args()
    input_matrix_path = args.input_matrix_path
    output_dir = args.output_dir
    delimiter = args.delimiter

    if delimiter == "t":
        expmatsep = "\t"
    else:
        expmatsep = ","

    Search_fields = ["sample_description", "dev_stage","study_title","sample_title", "tissue_lib", "tissue_type"]

    from data_processing import expression_matrix, metadata, read_write

    print(f"Getting SRR identifiers from expression matrix located at {input_matrix_path}...")
    sample_names = expression_matrix.load_only_samples(input_matrix_path, expmatsep = expmatsep)

    sub_outdir = os.path.join(output_dir, "Lemma_annotations")
    metadata_path = os.path.join(sub_outdir, "metadata.tsv")
    read_write.establish_dir(sub_outdir, isdir=True)
    print(f"Downloading metatadata of {len(sample_names)} RNA-seq samples to {metadata_path}...")
    metadata.download_metadata_for_exp_mat(sample_names, Search_fields, metadata_path)