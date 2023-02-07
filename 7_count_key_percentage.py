# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 09:25:33 2021

@author: titli

CODE FOR COUNTING PERCENTAGE OF EACH KEYS IN DATASET (HOW MANY PROTEIN IT BELONGS TO)
"""
import os, glob, time
#import multiprocessing
#from joblib import Parallel, delayed
import pandas as pd
import numpy as np 
import argparse

start_time = time.time()

parser = argparse.ArgumentParser(description='Key calculation - 1D to 3D.')
parser.add_argument('data_dir', type=str, help='Enter data dir')
parser.add_argument('filename', type=str, help='Enter the input filename')
parser.add_argument('total_protein', type=int, help='Enter total number of proteins in the dataset')
args = parser.parse_args()

data_dir = args.data_dir #'.\\..\\Dataset\\lexiographic\\mirror\\identified_mirrors\\d1d2\\0_Output_match_keys_flex0_plusMinus_1\\' 
total_protein = args.total_protein #43
filename = args.filename #'result_combined_all_triplets.csv' 
print(data_dir,"\t", filename, "\t", total_protein)

def count_percentage(filename, total_protein):
    df = pd.read_csv(data_dir+filename, header=0, sep=",")
    '''print("Original df len= ", len(df))
    keydict = {}
    for key in df['key'].tolist():
        if abs(key) not in keydict:
            keydict[abs(key)] = [key]
        else:
            keydict[abs(key)].append(key)
    keep_keys = []
    for k, v in keydict.items():
        v = list(set(v))
        if len(v) == 2:
            print(k, len(v), v)
            keep_keys.extend(v)
    #print(keep_keys)
    
    df = df[df['key'].isin(keep_keys)]
    print("keep df len=", len(df), "keep_keys=", len(keep_keys))
    '''
    df['key1'] = df['key'].abs()

    df['d1d2_sorted'] = ['_'.join(str(value) for value in sorted(tup)) for tup in zip(df['d1'], df['d2'])]
    df['triplets_sorted'] = ['_'.join(str(value) for value in sorted(tup)) for tup in zip(df['aa0'], df['aa1'], df['aa2'])]
    #print(df)
    
    df_grouped = df.groupby(['key1', 'd1d2_sorted', 'triplets_sorted'])#["protein"].apply(set).reset_index(name='protset')
    #print(df_grouped)
    
    with open (data_dir+'result_key_percentage.csv', 'w') as outfile:
        # write header
        outfile.writelines("key\td1d2_sorted\ttotal_prot\ttotal_unique_prot\tpercentage(total_unq/total_dataset)\n")
        for item in df_grouped: # get percentage of set of protein for every key,d1d2 combination 
            keycols = item[0]
            protcount = item[1]['protein'].tolist()
            protcountset = set(protcount)
            percentage = np.round((len(protcountset)/total_protein)*100)
            if percentage >= 80:
                print("\n", keycols[0], keycols[1], len(protcount), len(protcountset), percentage)
            outfile.writelines(str(keycols[0])+"\t"+str(keycols[1])+"\t"+str(len(protcount))+"\t"+ str(len(protcountset))+"\t"+str(percentage)+"\n")

# FUNC CALL
count_percentage(filename, total_protein)