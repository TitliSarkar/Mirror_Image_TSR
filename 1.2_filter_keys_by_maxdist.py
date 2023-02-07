# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 10:21:44 2021

@author: titli
"""
import os, glob, argparse
import math
import multiprocessing
from joblib import Parallel, delayed
import pandas as pd

#python code.py flexibility /path/to/input/files/    
parser = argparse.ArgumentParser()
parser.add_argument("maxdist_filter", type=int, default=15,help='Enter maxdist filter value') #the output will be saved in the same location
parser.add_argument("data_dir", type=str, default='./../Dataset/lexiographic/',help='Enter the choice of input folder location') #the output will be saved in the same location
parser.add_argument("extension", type=str, default='.triplets_29_35',help='Enter extension of the input files') #the output will be saved in the same location

#parser.add_argument("choiceofgrouping", type=str,help='Enter type of the operation from this list(lexiographic/aa_grp_0):')
#parser.add_argument("choiceofmirror", type=str,help='Enter type of the operation from this list(mirror/nomirror):')
#parser.add_argument("choiceoffilter", type=str,help='Enter type of the operation from this list(lexiographic/aa_grp_0):')
args = parser.parse_args()

# create directories
data_dir=args.data_dir #"./../Dataset/" 
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
subdir=data_dir+ "maxdist"+str(args.maxdist_filter)+"/"
if not os.path.exists(subdir):
    os.makedirs(subdir)
     
# read files
data_location = len(data_dir.split("/"))-1
print(data_location, "\n", data_dir.split("/"), "\n", subdir.split("/"))
files= [f.split("/")[data_location].split(".")[0] for f in glob.glob(data_dir+"*"+args.extension)]
print(len(files), files)

# filter keys
for file in files:
    df = pd.read_csv(data_dir+file+args.extension, sep="\t", header=0)
    df_filtered = df[df['maxDist'] <= args.maxdist_filter]
    df_filtered.to_csv(subdir+file+args.extension+'_maxdist'+str(args.maxdist_filter), index=False, sep="\t")
