#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 19:19:13 2019

@author: C00222141
"""
import os,glob
import pandas as pd
import argparse

df= pd.read_csv('sample_details.csv') # INPUT : keep this file at same location wtih the code
group_dict = df.groupby('group')['protein'].apply(list).to_dict()
for k,v in group_dict.items():
    print(k, len(v))
    
parser = argparse.ArgumentParser()
parser.add_argument("data_dir", help="enter the input files location", type=str)
args = parser.parse_args()

data_dir = args.data_dir # input files location


def extract_keys(group, files):
    file_count = len(files)
    key_dict = {}
    for f in files:
        df = pd.read_csv(data_dir+f+'.keys_29_35', sep="\t", usecols=[0], names=['key'])
        keyList = df['key'].tolist()
        print(len(keyList))
        for key in keyList:
            if key not in key_dict:
                key_dict[key] = [f]
            else:
                key_dict[key].append(f)
                
    keys_dict_count = {}         
    for k,v in key_dict.items():
        keys_dict_count[k]=len(v)
    print(keys_dict_count)
    df_out = pd.DataFrame(list(keys_dict_count.items()), columns=['key', '#files_present'])
    df_out['key_percent'] = round(df_out['#files_present']/file_count*100, 2)
    #print(df_out)
    
    df_only_100 = df_out[df_out['key_percent']==100]
    #print(df_only_100)
    
    df_out.to_csv(data_dir+'keys_percents_'+group+'_0.csv')
    df_only_100.to_csv(data_dir+'keys_percents_'+group+'_100.csv')
    

for group,files in group_dict.items():
    #if group != 'PKC':
        #continue
    print(group, len(files), files)
    extract_keys(group,files)
    


