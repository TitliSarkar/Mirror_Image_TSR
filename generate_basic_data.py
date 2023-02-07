# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 14:40:39 2019

@author: Titli
"""

import os, glob
import pandas as pd
import time
import numpy as np

data_dir = './../Dataset/'
subdir='lexiographic/'

filesList= sorted([f.split("/")[3].split(".")[0] for f in glob.glob(data_dir+"*.pdb")])
print(len(filesList))
for f in filesList:
    print (f)

df= pd.read_excel('Protein ID for TSR Manuscript 8-24-2018.xlsx',sheet_name='PKABC')
chain_dict = df.groupby('Chain')['PDB IDs'].apply(list).to_dict()
for k,v in chain_dict.items():
    print (k,len(v), v)

def countAminos(chain_dict, filesList):
    amino_dict ={}
    for chain,files in chain_dict.items():
        chainId = chain.upper()
        for fileName in files:
            inFile=open(data_dir+fileName+'.pdb','r')
            counter=0
            for i in inFile:
                if ((i[0:6].rstrip()=="ENDMDL") or (i[0:6].rstrip()=='TER' and i[21].rstrip()==chainId)):
                    break
                
                if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
                    break
                
                if(i[0:4].rstrip())=="ATOM"and(i[13:15].rstrip())=="CA"and(i[16]=='A'or i[16]==' ' and i[21:22]==chainId)and i[17:20]!= "UNK" :
                    counter+=1
                    
            amino_dict[fileName] = int(counter)
    
    amino_count = []
    for i in filesList:
        amino_count.append(int(amino_dict[i]))
        print(amino_dict[i])
    avg_aa = np.mean(amino_count)
    print("Average amino acids= ", avg_aa)

def countKeysTotal(filesList): # count total keys with freq
    print(" Function countKeysTotal():\n")
    total= []
    start=time.time()
    for f in filesList:
        inFile=open(data_dir+subdir+f+'.triplets_29_35','r') #commonTriplets_PKB_notOthers_29_35
        keys=0
        for i in inFile:
            keys +=1
        print(keys)
        total.append(int(keys))
    end=time.time()
    avg_total = np.mean(total)
    #print ("Average Total Keys (freq)= ", avg_total)
    print("Time taken= ", end-start)
    return avg_total

def countKeysUnique(filesList): # count total unique keys
    print(" Function countKeysUnique():\n")
    total= []
    start=time.time()
    for f in filesList:
        inFile=open(data_dir+subdir+f+'.keys_29_35','r') #commonTriplets_PKB_notOthers_29_35
        keys=0
        for i in inFile:
            keys +=1
        print(keys)
        total.append(int(keys))
    end=time.time()
    avg_total = np.mean(total)
    #print ("Average Total Keys (unq)= ", avg_total)
    print("Time taken= ", end-start)
    return avg_total

    
# func calls
countAminos(chain_dict, filesList)
print ("Average Total Keys (unq)= ", countKeysUnique(filesList))
print ("Average Total Keys (freq)= ", countKeysTotal(filesList))
