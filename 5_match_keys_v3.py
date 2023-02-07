# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 22:30:24 2019

@author: Titli
"""
import os
import pandas as pd
import numpy as np
import time
import argparse
import multiprocessing
from joblib import Parallel, delayed

start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("flexibility", default=0, help="enter the flexibility", type=int)
parser.add_argument("plusMinus", default=0, help="enter the plus/minus parameter", type=int)
parser.add_argument("path", default = './../Dataset/lexiographic/d1d2', help="enter path to the input file", type=str) #the output will be saved in a folder 'output_corr_seq_TSR' in the same location
args = parser.parse_args()

flexibility = args.flexibility 
plusMinus = args.plusMinus
print(" Flexibility = ", flexibility, " Plus/Minus= ", plusMinus)

data_dir = args.path # input files location
print(data_dir)
out_dir = data_dir+'0_Output_match_keys_flex'+str(flexibility)+'_plusMinus_'+str(plusMinus)+'/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    
triplets_df = pd.read_csv('all_amino_combinations.csv', sep='\t', names=['aa0','aa1','aa2'])
triplets_df[['aa0','aa1','aa2']] = [sorted(i) for i in triplets_df[['aa0','aa1','aa2']].to_numpy()]
triplets_df['combined'] = triplets_df['aa0']+'_'+triplets_df['aa1']+'_'+triplets_df['aa2']
triplets = triplets_df['combined'].tolist()
print(len(triplets)) # 1540

files = [] # each triplet -> one file for all protiens
for triplet in triplets:
    file = data_dir+triplet+'/combined.csv'
    if os.path.exists(file):
        files.append(file)
#files = [data_dir+triplet+'/combined.csv' for triplet in triplets] # each triplet -> one file for all protiens 
print(len(files))

def find_keys_by_condition_match(file):
    location = len(file.split("/"))-2
    triplet_name = file.split("/")[location]
    outfile = out_dir+triplet_name+'_flex'+str(flexibility)+'_plusOrMinus'+str(plusMinus)+'_d1crossd2_withProts.csv'
    #if os.path.exists(outfile):
        #return
    print(file, triplet_name)
    #result_dict = {}
    result_dflist = []
    start = time.time()
    print("\n\nprocessing file: ", file, triplet_name)
    
    # read current file   
    df= pd.read_csv(file, header=0)
    #print(df)
    #for each (d1,d2) pair in the current file, try to match it with all other files
    for index,row in df.iterrows():
        #distance = sorted([int(row['d1']), int(row['d2'])])
        d1 = int(row['d1']) # distance[0]
        d2 = int(row['d2']) #distance[1]
        k = int(row['key'])
        print(d1, type(d1), d2, type(d2), k, type(k))
        
        flex_d1 = round(d1*flexibility/100)
        flex_d2 = round(d2*flexibility/100)
            
        range_d1 = list(range(d1-flex_d1, d1+flex_d1+1))
        range_d2 = list(range(d2-flex_d2, d2+flex_d2+1))
        range_key = list(range(k-plusMinus, k+plusMinus+1))
        
        #print("\n",k, range_key, d1, d2, flex_d1, flex_d2, range_d1, range_d2, )
        
        # updation logic
        #logic1 = (df['d1'].isin(range_d1) & df['d2'].isin(range_d2))
        logic2 = (df['d1'].isin(range_d2) & df['d2'].isin(range_d1))
        logic3 = (df['key'].isin(range_key))
        logic4 = (df['d1'] != df ['d2'])

        # select row
        ####df_selected = df[(logic1 | logic2) & logic3] # select rows (original)
        df_selected = df[logic3 & logic4 & logic2] # select rows
        #if len(df_selected) != 0:
            #print(df_selected, "\n")
        result_dflist.append(df_selected)

    combined_csv = pd.concat([matched_df for matched_df in result_dflist])
    #combined_csv = combined_csv.drop(['protein'], axis=1)
    combined_csv.drop_duplicates(keep='first', inplace=True)
    combined_csv = combined_csv.sort_values(by=['key', 'd1', 'd2', 'protein'])
    print(len(combined_csv))
    if len(combined_csv) == 0:
        print("\nEmpty result! location = ", outfile)
    else:
        combined_csv.to_csv(outfile, index=False, sep = ",")
    print(file, "\n Processing time(sec): ", time.time()-start)

# func call 
num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)   
Parallel(n_jobs=num_cores, verbose=50)(delayed(find_keys_by_condition_match)(file)for file in files) # 1540 or less
#for file in files:
     #find_keys_by_condition_match(file)

print("CODE END\t TIME=(min)", np.round((time.time()-start_time)/60,2))
