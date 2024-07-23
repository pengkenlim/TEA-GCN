#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans , BisectingKMeans
from sklearn.metrics import silhouette_score
from data_processing import read_write

#User functions

def iterate_over_krange(data, k_list, k_cluster_assignment_dict_path , silhouette_coefficients_dict_path, mode="kmeans" ,randomstate=42):
    """Run kmeans clustering over a range of K"""
    kmeans_kwargs = {"init": "k-means++", "n_init": 5,"max_iter": 1000,"random_state": randomstate} #remove random
    silhouette_coefficients = []
    k_cluster_assignment_dict={}
    centroids_dict={}
    for k in k_list:
        if mode == "kmeans" or k < 3:
            kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        elif mode == "bisectingkmeans":
            kmeans = BisectingKMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(data)
        silhouette_coefficients.append(silhouette_score(data, kmeans.labels_))
        centroids_dict[k]= kmeans.cluster_centers_
        k_cluster_assignment_dict[k]= kmeans.labels_
        read_write.to_pickle(k_cluster_assignment_dict , k_cluster_assignment_dict_path)
        read_write.to_pickle( silhouette_coefficients, silhouette_coefficients_dict_path)
        print(f"K-means iteration at k={k} complete.SC:{silhouette_score(data, kmeans.labels_)}\n")

    return k_cluster_assignment_dict , silhouette_coefficients, centroids_dict

def select_k(silhouette_coefficients,k_cluster_assignment_dict):
    """"select Ks where silhoette coefficients peaks"""
    selected_list = []
    score = 0
    k_list = list(k_cluster_assignment_dict.keys())
    for k, sc in zip(k_list, silhouette_coefficients):
        if sc < score:
            selected_list.append(prev_k)
        else:
            prev_k = k
        score = sc
    selected_list.append(k)
    selected_list = sorted(list(set(selected_list)))
    return selected_list

def select_k_window_dep(silhouette_coefficients,k_cluster_assignment_dict, window_size = 10):
    """"select Ks with maximum sc every window"""
    selected_list = []
    k_list = list(k_cluster_assignment_dict.keys())
    window = []
    K_window =[]
    for k, sc in zip(k_list, silhouette_coefficients):
        if k % window_size ==0:
            selected_list.append(K_window[window.index(np.max(window))])
            window=[]
            K_window=[]
        else:
            window.append(sc)
            K_window.append(k)
    return selected_list

def select_k_window(silhouette_coefficients,k_cluster_assignment_dict, window_size = 10):
    """"select Ks with maximum sc every window"""
    selected_list = []
    k_list = list(k_cluster_assignment_dict.keys())
    max_k , min_k = np.max(k_list) , np.min(k_list)
    window = []
    window_min_K = min_k
    for window_max_K in range(min_k+window_size, max_k, window_size):
        for k in k_list:
            if k >= window_min_K and k < window_max_K:
                window.append(silhouette_coefficients[k_list.index(k)])
        selected_list.append(k_list[silhouette_coefficients.index(np.max(window))])
        window_min_K = window_max_K
        window = []
    return selected_list