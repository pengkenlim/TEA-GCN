#setting sys.path for importing modules
import pandas as pd
import requests
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

def metadata_json_from_ena(Accession, fields=["run_accession","fastq_aspera","fastq_ftp","fastq_bytes","library_layout","experiment_title","study_title", "sample_title","sample_description", "tissue_type", "tissue_lib","experimental_factor","local_environmental_context", "temperature", "host_growth_conditions","host_phenotype", "host_status"]):
    fields_string = "%2C".join(fields)
    url =f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22{Accession}%22&fields={fields_string}&format=json"
    response = requests.get(url)
    if response.status_code != 200:
        return "ERROR"
    return response.json()[0]

def download_metadata_for_exp_mat(Samples, Search_fields, outputpath):
    with open(outputpath, "w") as fout:
        fout.write("Accession"+ "\t"+ "\t".join(Search_fields) + "\n" )
        for idx , Accession in enumerate(Samples):
            try:
                response_json = metadata_json_from_ena(Accession, fields = Search_fields)
                contents = []
                for field in Search_fields:
                    contents.append(response_json[field])
                fout.write(Accession+ "\t"+ "\t".join(contents) + "\n" )

            except:
                try: # retry
                    response_json = metadata_json_from_ena(Accession, fields = Search_fields)
                    contents = []
                    for field in Search_fields:
                        contents.append(response_json[field])
                    fout.write(Accession+ "\t"+ "\t".join(contents) + "\n" )
                except:
                    fout.write(Accession+ "\t"+ "\t".join(["" for i in range(len(Search_fields))]) + "\n" )
            if idx % 100 == 0:
                print(f"Metadata for {idx} samples processed.")