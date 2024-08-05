
<img src="https://github.com/user-attachments/assets/682c7b18-fdc6-4353-bf35-28806b296484" alt="banner" width="600" />

 ## TEA-GCN: Two-Tier Ensemble Aggregation Gene Co-expression Network 


<img src="https://github.com/user-attachments/assets/f31ae18f-5846-49d7-b597-3234a7035ab2" alt="banner" width="800"/>

## What does this pipeline do?
This pipeline generates high-quality Gene Co-expression Networks (TEA-GCN ) that capture tissue/condition-specific co-expression.

For in-depth explanation and evaluation of the TEA-GCN methodology, please refer to our [preprint](https://doi.org/10.1101/2024.07.22.604713) :newspaper:

## Navigation
* [Generate TEA-GCN from your transcriptomic dataset](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#generate-tea-gcn-from-your-transcriptomic-dataset) :arrow_forward:
  * [Step 1. Setting up](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-filed#step-1-setting-up)
  * [Step 2. Generating partitions for your dataset](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-2-generating-partitions-for-your-dataset)
  * [Step 3. Building TEA-GCN](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-3-building-tea-gcn)
  * [Step 4. Post-processing TEA-GCN](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-4-post-processing-tea-gcn)
    
* [Gene Function Prediction using TEA-GCN](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#gene-function-prediction-using-tea-gcn) :dart:
  * [Step 1. Generating Co-expression Neighbourhoods of your genes-of-interest](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-1-generating-co-expression-neighbourhoods-of-your-genes-of-interest)
  * [Step 2. GSEA using Google Colab notebook](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-4-experimental-context-discovery-using-google-colab-notebook)
    
* [Discover experimental contexts underpinning TEA-GCN co-expression edges](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#discover-experimental-contexts-underpinning-tea-gcn-co-expression-edges) :mag_right: :zap::test_tube:
  * [Step 1. Downloading Metadata of RNA-seq samples from the European Nucleotide Archive (ENA)](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-1-downloading-metadata-of-rna-seq-samples-from-the-european-nucleotide-archive-ena)
  * [Step 2. Annotating Partitions with overrepresented lemmas](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-2-annotating-partitions-with-overrepresented-lemmas)
  * [Step 3. Generating Partition Rankings for your edges-of-interest](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-3-generating-partition-rankings-for-your-edges-of-interest)
  * [Step 4. Experimental context discovery using Google Colab notebook](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-4-experimental-context-discovery-using-google-colab-notebook)
    
* [Evaluating TEA-GCN Performance](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#evaluating-tea-gcn-performance) :chart_with_upwards_trend:
  * [Step 1. Preparing positive and negative edges](https://github.com/pengkenlim/TEA-GCN/tree/main?tab=readme-ov-file#step-1-preparing-positive-and-negative-edges)
  * [Step 2. Calculating ROC and PRC performance](https://github.com/pengkenlim/TEA-GCN/)


* [Citing TEA-GCN](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#citing-tea-gcn)
* [Complementary tools](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#complementary-tools)
* [Author and contact information](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#Author-and-contact-information)

## Generate TEA-GCN from your transcriptomic dataset

<img  src="https://github.com/user-attachments/assets/1e2d62ab-5a93-484e-b0da-6406fc9d1122" alt="banner" width="500"/>

### Step 1. Setting up
#### Clone repository to your local machine
```
$ git clone https://github.com/pengkenlim/TEA-GCN.git
```
#### Create an environment, and install packages
```
$ cd TEA-GCN
$ virtualenv -p python3 .venv
$ source ./.venv/bin/activate
$ pip install --upgrade pip
$ pip install -r ./setup/requirements.txt
```

#### Activate environment for subsequent use
```
$ cd TEA-GCN
$ source ./.venv/bin/activate
```
### Step 2. Generating partitions for your dataset
The TEA-GCN method uses of k-means clustering algorithm to divide gene expression data into partitions before gene co-expression determination. Expression data must be provided in the form of an expression matrix where 1) expression abundances are in the form of Transcript per Million (TPM), 2) rows correspond to different genes, and 3) columns correspond to different samples. The first row of the input expression matrix must consist of column headers (i.e. sample names) and the first column must consist of unique indices (i.e. gene identifiers).

#### Downloading sample data
We provide sample data of a gene expression matrices containing 500 and 5,000 _Arabidopsis thaliana_ public RNA-seq samples, respectively. (used in our [preprint](https://doi.org/10.1101/2024.07.22.604713) ).

You can download it from your browser using these links: 

1. [500 _Arabidopsis thaliana_ public RNA-seq samples](https://drive.google.com/file/d/1EuJPH772bsZaFqvQRzlcz7WptMCvOyjU/view?usp=sharing)
2. [5,000 _Arabidopsis thaliana_ public RNA-seq samples](https://drive.google.com/file/d/1EoLwtAoyc_xH2f_S-cyV8xt6E9sisuGa/view?usp=sharing)

#### Simplest implementation
```
$ python main/Generate_partitions.py --output_dir /path/to/output_directory --input_matrix_path /path/to/taxid3702_500n_expression_matrix.tsv
```

#### Full options
```
usage: Generate_partitions.py [-h] [-ks] [-ke] [-st] [-t] -o  -im  [-de] [-pc] [-rs]

Generate_partitions.py Load expression matrix. Perform sample standardization to correct batch effect. Embed samples from gene-space into PC-space. Perform k-means clustering across a
range of K. Outputs principal component embeddings of samples (PCA data), mean silhouette coefficients from each clustering iteration, k-means cluster assignments of samples in pickle
format.

options:
  -h, --help            show this help message and exit
  -ks , --k_start       Starting k for K-means clustering. The default is 10.
  -ke , --k_end         Ending k for K-means clustering. The Default is 50
  -st , --step          Step to increase between each k. The user is suggested to increase the step if the range is very big to decrease search space.
  -t , --threads        Number of threads to use for PCA and Kmeans clustering. Handled by Sklearn's parallelization. Try to specify lower than 64 for now.
  -o , --output_dir     Working directory to output data.
  -im , --input_matrix_path
                        Path of expression matrix to input
  -de , --delimiter     Delimiter for expression matrix. -de="t" for tab-separated (.tsv). -de="c" for comma separated (.csv). TSV by default.
  -pc , --pc_no         Number of Principal components to keep. Default = 100. It cannot be more than the number of samples and the number of genes in the expression matrix.
  -rs , --randomstate   randomstate (seed) for random seeding during k-means clustering. By default, 42 will be used.
```

### Step 3. Building TEA-GCN

After data partitioning, you can start building TEA-GCN by determining co-expression strength between every gene pair. Said co-expression strength is calculated based on measured correlation coefficients between genes from every dataset partition. If you would like more information, please refer to our [preprint](https://doi.org/10.1101/2024.07.22.604713).

#### Simplest implementation
```
$ python main/Run_TEA-GCN.py  --output_dir /path/to/output_directory --input_matrix_path /path/to/taxid3702_500n_expression_matrix.tsv
```
Intermediate TEA-GCN will be output to `/path/to/output_directory/RAW_GCN/<METHOD>/<NETWORK_NAME>`. `<NETWORK_NAME>` is generated automatically based on the `<METHOD>_<K>k_<AGGREGATION_FUNCTION>` format (e.g. TEA_8k_RAvg)

#### Full options
```
usage: Run_TEA-GCN.py [-h] [-w] [-t] -o  [-de] [-cc] [-am] [-k] -im

Run_TEA-GCN.py. Construct TEA-GCN from expression data (in the form of an expression matrix). Requires prior data partitioning using Generate_partitions.py in order to run.

options:
  -h, --help            show this help message and exit
  -w , --workers        Number of workers for parallelization. Affects precalc of all coefficients. But only affects the calc of bicor.
  -t , --threads        Number of threads for numpy linear algebra operations. Affects PCC and SCC calc.
  -o , --output_dir     Directory to output. Must be the same as for Generate_partitions.py
  -de , --delimiter     Delimiter for expression matrix. -de="t" for tab-separated (.tsv). -de="c" for comma-separated (.csv). TSV by default.
  -cc , --correlation_coefficient
                        select 'TEA' to enable Coefficient Aggregation. Select "PCC","SCC","bicor" for constructing Partition Aggregation-only GCN
  -am , --aggregation_method
                        The default is RAvg.
  -k , --k_clusters     Number of clusters to partition the expression data. For k=0, will use best k as determined in Generate_partitions.py step
  -im , --input_matrix_path
                        Path of expression matrix to input
```

### Step 4. Post-processing TEA-GCN

After building a TEA-GCN that describes the co-expression strengths between every gene pair (edges). We can then load back in all edges for rank transformation (into Mutual Ranks) and generate standardized ranks for meaningful cross-comparison of TEA-GCN across different species.

#### Simplest implementation

```
$ python ./main/Rank_transform.py --output_dir /path/to/output_directory
```

Completed TEA-GCN will be located in `/path/to/output_directory/Completed_GCN/<NETWORK_NAME>`.

#### Full options

```
usage: Rank_transform.py [-h] -o  [-n] [-d] [-full]

Rank_transform.py calculates Mutual Ranks (MR) and Highest reciprocal ranks (HRR) and their standardized forms raw correlation scores. Rank_transform.py outputs the final complete TEA-GCN
as directory separate from the intermediate TEA-GCN generated by Run_TEA_GCN.py

options:
  -h, --help            show this help message and exit
  -o , --output_dir     Working directory containing required data of which to output data. Must be same as output_dir for Run_TEA_GCN.py and Generate_partitions.py
  -n , --net_name       name of the network to calculate MR and HRR for. By default, all networks will be calculated
  -d , --delete_net     If set to True, Rank_transform.py will delete intermediate TEA-GCN generated by Run_TEA_GCN.py. False by default.
  -full , --full_attributes
                        If set to True, Rank_transform.py will output a full set of attributes for every edge which includes: raw co-exp. strengths, MR, HRR, and their Z score variants.
                        False by default where only raw co-exp, MR, and Z(MR) are reported.
```

#### Data structure of Network

The resultant network is contained within a directory consisting of files that each describe edges connecting to each gene (source gene). The files are named after genes. 

Example of one file:

```
$ less /path/to/<NETWORK_NAME>/AT5G65360.1
```

<img src="https://github.com/user-attachments/assets/bde075ec-7015-4661-9caf-6d477cb50c94" alt="banner"  width="600"/>

Description of columns:

* `Target`
 * Describes the target gene of the edge 
* `Co-exp_Str`
 * Co-expression Strength of edge calculated by `Run_TEA-GCN.py`
* `Co-exp_Str_MR`
  * Co-expression Strength of edge transformed into Mutual Ranks (MR)
* `zScore(Co-exp_Str_MR)`
  * z-score standardized MR -- Helpful for cross-comparing GCNs

Back to [Navigation](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#navigation)

## Gene Function Prediction using TEA-GCN

<img src="https://github.com/user-attachments/assets/af7a133e-bad2-4a72-852f-add13fb173bb" alt="banner"  width="700"/>

### Step 1. Generating Co-expression Neighbourhoods of your genes-of-interest

#### Simplest implementation

```
$ python ./main/Prepare_neighbourhoods.py --output_dir /path/to/output_directory --gene_list AT3G02480.1
```

#### Full Options

```
usage: Prepare_neighbourhoods.py [-h] -o  -g  [-n]

Prepare_neighbourhoods.py. Prepare co-expression neighbourhoods of Genes-of-intrests from completed TEA-GCNs for gene function prediction using the GSEA algorithm.

options:
  -h, --help          show this help message and exit
  -o , --output_dir   Directory to output. Must be the same as for Generate_partitions.py, Run_TEA-GCN.py, and Rank_transform.py.
  -g , --gene_list    list of comma-separated genes. spaces will be removed. E.g. 'GENE1,GENE2,GENE3'
  -n , --net_name     name of the TEA-GCN to extract. Set to "auto" if you only generated one TEA-GCN. "auto" by default
```

#### Format of Co-expression Neighbourhood files

Example of one file:

```
$ less /path/to/output_directory/Co-exp_Neighbourhoods/AT5G65360.1
```

<img src="https://github.com/user-attachments/assets/31c0d11a-ccbc-4d8a-8f18-acad1c7e4899" alt="banner"  width="300"/>

Description of columns:

* `Target`
 * Describes the target gene of the edge 
* `Score`
 * z-score standardized MR determined by `Rank_transform.py`
* `Centered_rank`
  * `Score` ranked inversely (higher values, have larger ranks) and centered (median rank = 0).

### Step 2. GSEA using Google Colab notebook

Refer to the Colab notebook (linked below) to predict the biological function of your gene-of-interest based on their co-expression neighbourhood:

[TEA-GCN_Function_pred_using_GSEA.ipynb](https://colab.research.google.com/drive/1Qjog6kc5QtVOl0Rdh8fhnZDVBF4sCz8h?usp=sharing)

![note_1](https://github.com/user-attachments/assets/e9e1d64b-224a-4bba-b0f4-08f4ea978264)

![note_2](https://github.com/user-attachments/assets/3aef94db-8da5-4c24-a921-084c67971357)

![note_3](https://github.com/user-attachments/assets/c42e241e-23d2-4b76-93f3-53a5d2c9414e)

![note_4](https://github.com/user-attachments/assets/a1223e67-ac0a-4f37-ab9f-687528780c3a)

Back to [Navigation](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#navigation)

## Discover experimental contexts underpinning TEA-GCN co-expression edges

<img  src="https://github.com/user-attachments/assets/bf40485c-768d-4eb7-be0f-4ed39530b884" alt="banner" width="700"/>

### Step 1. Downloading Metadata of RNA-seq samples from the European Nucleotide Archive (ENA)

This script downloads the metadata of RNA-seq samples from the European Nucleotide Archive (ENA). 

The metadata downloaded includes these search fields*:

* `sample_description`
* `dev_stage`
* `study_title`
* `sample_title`
* `tissue_lib`
* `tissue_type`

_*for more information regarding search fields and metadata hosted by ENA, click_ [_here_](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/metadata.html)

#### Simplest implementation

```
$ python ./main/Download_metadata.py --output_dir /path/to/output_directory --input_matrix_path /path/to/taxid3702_500n_expression_matrix.tsv
```

#### Full options

```
usage: Download_metadata.py [-h] -o  [-de] -im

Download_metadata.py Download metadata of SRR accessions from European Nucleotide Archive (ENA).

options:
  -h, --help            show this help message and exit
  -o , --output_dir     Directory to output. Must be the same as for Generate_partitions.py, Run_TEA-GCN.py and Rank_transform.py.
  -de , --delimiter     Delimiter for expression matrix. -de="t" for tab-separated (.tsv). -de="c" for comma separated (.csv). TSV by default.
  -im , --input_matrix_path
                        Path of expression matrix to input
```

### Step 2. Annotating Partitions with overrepresented lemmas

<img  src="https://github.com/user-attachments/assets/ed3f63e4-c69a-45ea-90da-9dfc1419180c" alt="banner" width="700"/>

#### Simplest implementation

Download and install your spaCy model of choice for Natural Language Processing (NLP). Below we chose the `en_core_web_sm` model. This model will be used to lemmatize the metadata of RNA-seq samples downloaded from ENA.

```
$ python -m spacy download en_core_web_sm
```

Running this script will find overrepresented lemmas to annotate partitions

```
$ python ./main/Annotate_partitions.py  --output_dir /path/to/output_directory --input_matrix_path /path/to/taxid3702_500n_expression_matrix.tsv
```

#### Full options

```
usage: Annotate_partitions.py [-h] -o  [-de] -im  -k  [-m] [-p]

Annotate_partitions.py Lemmatize metadata of RNA-seq samples. Annotate dataset partitions with overrepresented lemmas

options:
  -h, --help            show this help message and exit
  -o , --output_dir     Directory to output. Must be the same as for Generate_partitions.py, Run_TEA-GCN.py, and Rank_transform.py.
  -de , --delimiter     Delimiter for expression matrix. -de="t" for tab-separated (.tsv). -de="c" for comma separated (.csv). TSV by default.
  -im , --input_matrix_path
                        Path of expression matrix to input
  -k , --k_clusters     Number of clusters to partition the expression data. For k=0, will use best k as determined in Generate_partitions.py step.
  -m , --model          spaCy model to use for lemmatization. 'en_core_web_sm' will be used by default.
  -p , --permutations   Number of permutations to shuffle during the calculation of overrepresentation p-values. Default = 10,000

```

### Step 3. Generating Partition Rankings for your edges-of-interest

<img  src="https://github.com/user-attachments/assets/9252f188-670b-4d95-8a99-0d86dc19527d" alt="banner" width="250" />

#### Defining edges-of-interest

Prepare edges-of-interest as input like so.

```
$ less /path/to/edges_of_interest.csv
```

<img  src="https://github.com/user-attachments/assets/a4d1fa68-765c-42b1-96a0-4f2c2c553e20" alt="banner" width="210" />

Each row describes one edge. Two comma-separated columns for the two genes that connect the edge.

#### Simplest implementation

```
python ./main/Get_Partition_Rankings.py -k 93 --input_matrix_path /path/to/taxid3702_500n_expression_matrix.tsv --edges_path /path/to/edges_of_interest.tsv --output_dir /path/to/output_directory
```

#### Full options

```
usage: Get_Partition_Rankings.py [-h] [-w] [-t] -o  [-de] [-k] -im  -ep  [-on]

Get_Partition_Rankings.py.

options:
  -h, --help            show this help message and exit
  -w , --workers        Number of workers for parallelization. Affects precalc of all coefficients. But only affects the calc of bicor.
  -t , --threads        Number of threads for numpy linear algebra operations.
  -o , --output_dir     Directory to output.
  -de , --delimiter     Delimiter for expression matrix. -de="t" for tab-separated (.tsv). -de="c" for comma separated (.csv). TSV by default.
  -k , --k_clusters     Number of clusters to partition the expression data. For k=0, will use best k as determined in Generate_partitions.py step.
  -im , --input_matrix_path
                        Path of expression matrix to input
  -ep , --edges_path    Path to list of edges to calculate.
  -on , --outname_prefix
                        prefix of to output calculations.
```

### Step 4. Experimental context discovery using Google Colab notebook

<img  src="https://github.com/user-attachments/assets/1357907f-a701-4585-b522-9223c4c2d14b" alt="banner" width="700"/>

#### Removing irrelevant lemmas

To discover experimental contexts underpinning gene co-expression, we use a GSEA-like approach to detect overrepresented lemmas that annotate partitions with high co-expression. However, not all overrepresented lemmas are informative or relevant to experimental conditions. Including them in the downstream analysis will only confound the interpretation of experimental contexts and decrease the sensitivity via harsher multiple-testing correction.

[Step 2](https://github.com/pengkenlim/TEA-GCN) generates a list of overrepresented lemmas that are detected in your dataset partitions in `/path/to/output_directory/Lemma_annotations/overrepresented_lemmas.tsv`. You can remove irrelevant lemmas from this file. It will be used as input in the Colab notebook below to exclude them from the analysis.

```
$ less /path/to/output_directory/Lemma_annotations/overrepresented_lemmas.tsv
```

<img  src="https://github.com/user-attachments/assets/874aa637-1dfd-49af-91cf-f70ab597e004" alt="banner" width="100"/>


#### Colab notebook

Refer to the Colab notebook (linked below) to uncover the Experimental contexts underpinning gene co-expression:

[Uncover_the_Experimental_contexts_Underpinning_Gene_Co-expression.ipynb](https://colab.research.google.com/drive/1DKVj8lNCnYXbhseMeIJzQk9TBGy1lNy-?usp=sharing)

<img  src="https://github.com/user-attachments/assets/ad84e67f-e8cd-4e37-81c2-bf05b8478fee" alt="banner" width="700"/>

Back to [Navigation](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#navigation)

## Evaluating TEA-GCN Performance

### Step 1. Preparing positive and negative edges

**COMING SOON**

Back to [Navigation](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#navigation)

## Citing TEA-GCN

You can cite our preprint for now :)

```
Lim, P. K., Wang, R., Velankanni, J. P. A., & Mutwil, M. (2024).
Constructing Ensemble Gene Functional Networks Capturing Tissue/condition-specific Co-expression from Unlabled Transcriptomic Data with TEA-GCN
(p. 2024.07.22.604713).
bioRxiv. https://doi.org/10.1101/2024.07.22.604713
```
## Complementary tools

**COMING SOON**

Back to [Navigation](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#navigation)

## Author and contact information


### Marek Mutwil :bearded_person:

<img  src="https://github.com/user-attachments/assets/21bc4e2c-ee6d-4534-a69b-b9d735322df7" alt="banner" width="100"/>\
_I am the project supervisor. I advise the students and buy them meals sometimes._

Principal Investigator and Assoc. Prof.,\
[School of Biological Sciences, Nanyang Technological University](https://www.ntu.edu.sg/sbs)\
mutwil@ntu.edu.sg:email:\
[Lab website](https://www.plant.tools/):link:


### Peng Ken Lim :boy:

<img  src="https://github.com/user-attachments/assets/439a32b9-7e34-4a65-98e4-bc9616956ef9" alt="banner" width="120"/>\
_I am the project co-supervisor. I wrote most of the code._

PhD student,\
[School of Biological Sciences, Nanyang Technological University](https://www.ntu.edu.sg/sbs)\
pengken001@e.ntu.edu.sg:email:

### Ruoxi Wang :baby:

_I wrote some of the code and conducted comparative-GCN analysis for this project._

Undergraduate student,\
[Shanghai Jiao Tong University](https://en.sjtu.edu.cn/)

### Jenet Princy Antony Velankanni :baby:

_I helped with mining Human and Yeast gene annotations for this project._

Undergraduate student,\
[School of Biological Sciences, Nanyang Technological University](https://www.ntu.edu.sg/sbs)

Back to [Navigation](https://github.com/pengkenlim/TEA-GCN/?tab=readme-ov-file#navigation)
