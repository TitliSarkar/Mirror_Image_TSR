#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:51:10 2019

@author: c00222141
"""
import glob, os, time
import pandas as pd 
import argparse
import multiprocessing
from joblib import Parallel, delayed
from itertools import chain
import numpy as np

start_time = time.time()

#python code.py path/to/input/files/    
parser = argparse.ArgumentParser()
parser.add_argument("data_dir", help="enter path to the input file", type=str) 
args = parser.parse_args()


data_dir = args.data_dir # input files location
print(data_dir)
data_location = len(data_dir.split("/"))-1
print(data_dir.split("/"))
print(data_location)

files = [x.split("/")[data_location].split(".")[0].split("_")[0] for x in glob.glob(data_dir+"*_distance_added.csv")]
print(len(files), files)

def seperate_triplets(file):
    print(file)
    #df= pd.read_csv(data_dir+file+"_distance_added.csv", sep="\t", usecols =[0,1,2,3,4,5,6,7,8,9],names=['key','amino0','pos0','amino1','pos1','amino2','pos2', 'd1','d2','protein'])
    df= pd.read_csv(data_dir+file+"_distance_added.csv", sep="\t", header=0, usecols =[0,1,2,3,4,5,6,7,8,9])

    df[['aa0','aa1','aa2']] = [sorted(i) for i in df[['aa0','aa1','aa2']].values]
    #print(df)
    df_grouped = df.groupby(['aa0', 'aa1','aa2']) 
    
    for name, group in df_grouped: 
        #print(name) 
        #print(group) 
        triplet_name = name[0]+"_"+name[1]+"_"+name[2]
        out_dir = data_dir+triplet_name+'/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        group.to_csv(out_dir+triplet_name+"_"+file+".csv", sep="\t", index=False)
        
# func call  
num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)
Parallel(n_jobs=num_cores, verbose=50)(delayed(seperate_triplets)(fileName)for fileName in files)
#for file in files:
    #seperate_triplets(file)
print("CODE END\t TIME=(min)", np.round((time.time()-start_time)/60,2))
