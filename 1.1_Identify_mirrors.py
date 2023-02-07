#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 10:56:11 2019

@author: c00222141
"""
import os, glob
import argparse
import pandas as pd 
from joblib import Parallel, delayed
import multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument("data_dir", default='./../Dataset/lexiographic/mirror/', help="enter path to the input file", type=str) # ./../Dataset/Intra_data_negetive_key/
parser.add_argument("extension", default='.triplets_29_35_mirror', help="enter file extension") 
args = parser.parse_args()

data_dir = args.data_dir #
out_dir = data_dir+'identified_mirrors/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    
data_location = len(data_dir.split("/"))-1
print(data_dir.split("/"))
print(data_location)

files = [x.split("/")[data_location].split(".")[0] for x in glob.glob(data_dir+"*"+args.extension)]
print(len(files), files)

def identify_mirrors_from_negetive_keys(file):
    print(file)
    #df = pd.read_csv(data_dir+file+'.triplets_29_35_aa_grp_0_mirror', header=0, sep="\t", usecols=[0,1,2,3,4,5,6], names=['key','amino0','pos0','amino1','pos1','amino2','pos2'])
    df = pd.read_csv(data_dir+file+args.extension, header=0, sep="\t")

    keys = list(set(df['key'].tolist())) # take set
    #keys = [-1,-1,1,1,1,1,1,4,5,6,-5,-5,-3,-5,3,2,1]
    #print(len(keys))
    
    seen = [] # list of keys having both -ve and +ve values: mirror candidates
    count_dict = {}
  
    #For each element of array 
    for i in keys:  
        if abs(i) not in count_dict:
            count_dict[abs(i)] = 1
        else: 
            seen.append(abs(i))
            count_dict[abs(i)] = 0
            
    #print(count_dict)      
    if len(seen) == 0:
        print("No mirros seen for protein: ", file)
    seen = sorted(seen)
    #print(seen)
    #print("#of mirror unq keys (one occurence per -ve/+ve) = ", len(seen))
    
    seen_unq = list(set(seen))
    seen_neg = [(-1)*abs(x) for x in seen_unq]
    seen_unq.extend(seen_neg)
    #print(seen_unq)
    
    df_seen = df[df['key'].isin(seen_unq)]
    temp = df_seen['key'].abs()
    df_seen['key_abs'] = temp
    #df_seen.assign(key_abs=df_seen['key'].abs())
    df_seen_sorted = df_seen.sort_values('key_abs', axis=0, kind='quicksort')
    df_seen_sorted1 = df_seen_sorted.drop(columns = ['key_abs'])
    df_seen_sorted1.to_csv(out_dir+file+args.extension+'_identified', index=False, sep='\t')
    
    '''# find mirrors of each key
    df['key_abs'] = df['key'].abs()
    grouped_df = df.groupby("key_abs") # dataframe group by "key"
    for g in grouped_df.groups.keys(): # process for each key
        print(grouped_df.get_group(g))
    #DataFrame.copy(deep=True)[source]'''

#func call
num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)
Parallel(n_jobs=num_cores, verbose=50)(delayed(identify_mirrors_from_negetive_keys)(fileName)for fileName in files)
#for file in files:
    #identify_mirrors_from_negetive_keys(file)
