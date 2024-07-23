<img src="https://github.com/user-attachments/assets/682c7b18-fdc6-4353-bf35-28806b296484" alt="banner" width="600"/>

# The TEA-GCN pipeline
## TEA-GCN: Two-Tier Ensemble Aggregation Gene Co-expression Network 

<img src="https://github.com/user-attachments/assets/f31ae18f-5846-49d7-b597-3234a7035ab2" alt="banner" width="800"/>

## What does this pipeline do?
This pipeline generates high-quality Gene Co-expression Networks (TEA-GCN ) that capture tissue/condition-specific co-expression.

## Navigation
* [Generate TEA-GCN from your transcriptomic dataset](https://github.com/pengkenlim/TEA-GCN/edit/main/README.md#generate-tea-gcn-from-your-transcriptomic-dataset)
  * [Step 1. Setting up](https://github.com/pengkenlim/TEA-GCN/edit/main/README.md#step-1-setting-up)
  * [Step 2. Generating partitions for your dataset](https://github.com/pengkenlim/TEA-GCN/edit/main/README.md#step-2-generating-partitions-for-your-dataset)
  * [Step 3. Building TEA-GCN](https://github.com/pengkenlim/TEA-GCN/edit/main/README.md#step-3-building-tea-gcn)
  * [Step 3. Post-processing TEA-GCN](www.google.com)
*[Gene Function Prediction using TEA-GCN](www.google.com)
  *[Step 1. Generating Co-expression Neighbourhood](www.google.com)
  *[Step 2. Prepare functional annotation gene sets](www.google.com)
  *[Step 2. GSEA using Google colab notebook](www.google.com)
*[Discover experimental contexts underpinning TEA-GCN co-expression edges](www.google.com)
  *[Step 1. Generating Partition Rankings for edges-of-interest](www.google.com)
  *[Step 2. Annotating Partitions with overrepresented lemmas](www.google.com)
  *[Step 3. Experimental context discovery using Google colab notebook](www.google.com)
* [Evaluating TEA-GCN Performance](www.google.com)
  *[Step 1. Preparing positive and negative edges](www.google.com)
  *[Step 2. Calculating ROC and PRC performance](www.google.com)

## Generate TEA-GCN from your transcriptomic dataset
### Step 1. Setting up
#### Clone repository to local machine
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
We provide sample data of a gene expression matrix containing 500 _Arabidopsis thaliana_ public RNA-seq samples (generated in **INSERT DOI**).

You can download it from your browser using this link: https://drive.google.com/file/d/1E0eJd6AsJw6VvXUOfDXuhUOjpgL0zk4u/view?usp=sharing

or via the command line:

```
$ wget -O /path/to/taxid3702_500n_expression_matrix.tsv https://drive.google.com/uc?id=1E0eJd6AsJw6VvXUOfDXuhUOjpgL0zk4u
```

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
After data partitioning, you can start building TEA-GCN by determining co-expression strength between every gene pair. Said co-expression strength is calculated based on measured correlation coefficients between genes from every dataset partition. For more information, refer to **INSERT DOI**.

#### Simplest implementation
```
$ python main/Run_TEA-GCN.py  --output_dir /path/to/output_directory --input_matrix_path /path/to/taxid3702_500n_expression_matrix.tsv
```

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
