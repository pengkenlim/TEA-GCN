#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         print(parent_module)
         sys.path.insert(0, parent_module)

from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np

def standardize_transform(Matrix, n_pcs=100):
    """normalize Matrix, transpose into dataframe"""
    #normalize colummns (within samples)
    #print("normalize colummns (within samples)")
    sample_normalized_matrix = pd.DataFrame(StandardScaler().fit_transform(Matrix).T,
                                            index=Matrix.columns,
                                            columns=Matrix.index)
    #print("PCA transformation")
    #PCA transformation
    #pca = PCA()
    pca = IncrementalPCA(n_components=n_pcs, batch_size=5000) # Vital to prevent memory overload
    pca_data=pca.fit_transform(sample_normalized_matrix)
    pca_data = pd.DataFrame(pca_data, index= Matrix.columns)
    #calculate total explained variance in the pcs that were kept
    pc_variances = np.round(pca.explained_variance_ratio_*100, decimals=1)

    return pca_data , pc_variances

def subset_pca_data(pca_data, n_pcs=1000):
    """subset pca_data to a set number of pcs"""
    pca_data = pca_data.iloc[:, 0:n_pcs ]
    return pca_data 