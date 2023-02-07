# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 09:47:04 2021

@author: titli

## COMBILE RESULTS FOR ALL TRIPLETE TO ONE FILE
"""
import os, glob, time
import pandas as pd
import time
import argparse
import multiprocessing
from joblib import Parallel, delayed
import numpy as np 

start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("data_dir", type=str, default='./../Dataset/lexiographic/', help="enter path to the input/output files")
args = parser.parse_args()
print(args.data_dir)

def combine_result_files(location):
    #triplets_df = pd.read_csv('all_amino_combinations.csv', sep='\t', names=['aa0','aa1','aa2'])
    #triplets_df[['aa0','aa1','aa2']] = [sorted(i) for i in triplets_df[['aa0','aa1','aa2']].to_numpy()]
    #triplets_df['combined'] = triplets_df['aa0']+'_'+triplets_df['aa1']+'_'+triplets_df['aa2']
    #triplets = triplets_df['combined'].tolist()
    #print(len(triplets)) # 1540

    combined_csv = pd.concat([pd.read_csv(f) for f in glob.glob(args.data_dir+"*.csv")]) ## combile results for all triplets in one file 
    #print(combined_csv)
    combined_csv.to_csv(args.data_dir+"result_combined_all_triplets.csv", index=False, encoding='utf-8-sig')

# func call
combine_result_files(args.data_dir)

print("CODE END. TIME TAKEN(min)=", np.round((time.time()-start_time)/60,2))