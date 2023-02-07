#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:51:10 2019

@author: c00222141

### Extract all keys for each amino acid triplet combinations
"""
import glob, os
import pandas as pd 
import argparse
import multiprocessing
from joblib import Parallel, delayed

#python code.py path/to/input/files/    
parser = argparse.ArgumentParser()
parser.add_argument("path", help="enter path to the input file", type=str) 
args = parser.parse_args()

data_dir = args.path # input files location
print(data_dir)

def combine_data(triplet):
    #print(triplet)
    all_filenames = [x for x in glob.glob(data_dir+triplet+"/"+triplet+"*.csv")] 
    #print(triplet, len(all_filenames))
    if len(all_filenames)!=0:
        #print(all_filenames[0])
        #combined_csv = pd.concat([pd.read_csv(f,sep="\t",names=['key','amino0','pos0','amino1','pos1','amino2','pos2', 'd1','d2','protein']) for f in all_filenames])
        #combined_csv = pd.concat([pd.read_csv(f,sep="\t",header=0) for f in all_filenames])
        li = []
        for f in all_filenames:
            if not os.path.exists(f):
                continue
            df = pd.read_csv(f,sep="\t",header=0)
            if len(df) > 0:
                li.append(df)
        print(triplet, len(all_filenames), len(li))
        if len(li) > 0:
            combined_csv = pd.concat(li, ignore_index=True)
            combined_csv.to_csv(data_dir+triplet+"/"+"combined.csv", index=False)
        
# func call  
triplets_df = pd.read_csv('all_amino_combinations.csv', sep='\t', names=['aa0','aa1','aa2'])
triplets_df[['aa0','aa1','aa2']] = [sorted(i) for i in triplets_df[['aa0','aa1','aa2']].values]
triplets_df['combined'] = triplets_df['aa0']+'_'+triplets_df['aa1']+'_'+triplets_df['aa2']
triplets = triplets_df['combined'].tolist()

num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)   
#Parallel(n_jobs=num_cores, verbose=50)(delayed(combine_data)(t)for t in triplets)
for t in triplets:
    combine_data(t)
print("Code end.")
