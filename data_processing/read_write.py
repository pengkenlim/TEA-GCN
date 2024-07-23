#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import pickle

# func definitions
def load_pickle(read_path):
    with open(read_path, 'rb') as f:
        data = pickle.load(f)
    return data

def to_pickle( data , write_path):
    with open(write_path, 'wb') as f:
        pickle.dump(data, f)

def establish_dir(path, isdir = False):
    """"create all directiories within the path if it does not already exist"""
    if isdir:
        path_dir = path
    else:     
        path_dir = "/".join(path.split("/")[:-1])
    if not os.path.exists(path_dir):
        os.makedirs(path_dir)
