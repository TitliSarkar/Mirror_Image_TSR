0# -*- coding: utf-8 -*-

import pandas as pd
import multiprocessing
from joblib import Parallel, delayed
import os, glob

data_dir = './../Dataset/lexiographic/'

filesList = sorted([x.split("/")[4].split(".")[0] for x in glob.glob(data_dir+"*.keys_29_35")])
print(len(filesList), filesList)

# Step 1
#finds unique common keys among given set of proteins
def find_common_keys(inputProts):
    keys_dict = {} # key=ProtKey, val=Protname

    for p in inputProts: 
        with open(data_dir+p+'.keys_29_35', 'r') as file: # .keys_29_35
            for r in file: # row is list
                if not r or len(r.strip().split()) != 2:
                    continue
                row = r.strip().split()
                #print (row[0], row[1])
                if row[0] not in keys_dict.keys():   # row[0] = prot, row[1] = key
                    keys_dict[row[0]] = []
                    keys_dict[row[0]].append(p)
                else:
                    keys_dict[row[0]].append(p)

    # remove duplicate proteins coming from multiple instances of a key within same
    for key_, val_ in keys_dict.items():
        val_ = list(set(val_))
        keys_dict[key_]=val_

    # remove entries with not common keys
    keys_dict = { k:v for k,v in keys_dict.items() if (len(v)==len(inputProts) and sorted(v)==sorted(inputProts))} 
    
    # find common keys
    common_keys = list(keys_dict.keys())
    print(" #of unq common keys= ",len(common_keys))
    
    #write keys in a an output file
    with open(data_dir+"commonKeys.txt", "w")as f:
        for k,v in keys_dict.items():
            if(len(v) == len(inputProts)):
                f.write(str(k)+"\n")
                #print(k,len(v))
    return common_keys

# Step 2: find common key details
def parallel_code(p, commonKeys):
    print(p)
    result_file = open(data_dir+p+'.commonTriplets_29_35','w') #.commonTriplets_29_35
    result_file1 = open(data_dir+p+'.commonKeys_29_35','w') #.commonKeys_29_35

    with open(data_dir+p+'.triplets_29_35', 'r') as file: # .triplets_29_35
            for row in file:
                if not row or len(row.strip().split()) != 20:
                    continue
                if row.strip().split()[0] in commonKeys:
                    result_file.writelines(row)

    with open(data_dir+p+'.keys_29_35', 'r') as file1: # .keys_29_35
            for row1 in file1:
                if row1.split()[0] in commonKeys:
                    result_file1.writelines(row1)

def find_commonkey_details(inputProts, commonKeys):
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores, verbose=50)(delayed(parallel_code)(fileName,commonKeys)for fileName in inputProts)

# func call    
find_common_keys(filesList)
commonKeys = []
fc = open(data_dir+"commonKeys.txt", "r")
for line in fc:
    if not line:
        continue
    commonKeys.append(line.strip())
print(commonKeys)

find_commonkey_details(filesList, commonKeys)
