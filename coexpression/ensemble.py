#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import numpy as np
import scipy.stats
import warnings

def nan2zeros(array):
    array_2 = array.copy()
    array_2[np.isnan(array_2)] = 0
    return array_2

def rectify(array):
    array_2 = array.copy()
    array_2[array_2 < 0] = 0
    return array_2

def scale(cor):
    return (cor +1)/2

def unscale(cor):
    return (cor*2)-1

def get_ranks(array,axis):
    rev_ranks = (scipy.stats.rankdata(array, method="max", nan_policy= "omit", axis=axis)).astype(np.float64) #bigger the values have larger (lower) rank. Tied values will be assigned the max rank
    ranks =  rev_ranks*(-1) + np.nanmax(rev_ranks) +1
    return rev_ranks, ranks

def RRWA(array , ranks, axis):
    custom_weights = 1/ranks
    return unscale(np.average(scale(array), weights = custom_weights, axis = axis))

def RWA(array , rev_ranks, axis):
    custom_weights = rev_ranks
    return unscale(np.average(scale(array), weights = custom_weights, axis = axis))

def aggregate(cor_values, mode, axis = 0):
    warnings.filterwarnings(action='ignore', message='Mean of empty slice')
    cor_values = np.array(cor_values)
    cor_values = nan2zeros(cor_values)
    rectified_cor_values = rectify(cor_values)
    rev_ranks, ranks = get_ranks(rectified_cor_values, axis)
    if mode == "RRWA":
        return RRWA(rectified_cor_values , ranks, axis)
    elif mode == "RWA":
        return RWA(rectified_cor_values, rev_ranks, axis)
    elif mode == "RAvg":
        return np.nanmean(rectified_cor_values, axis = axis)
    elif mode == "Avg" :
        return np.nanmean(cor_values, axis= axis)
    elif mode == "Max" :
        return np.nanmax(cor_values, axis = axis)