# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 18:41:06 2019

@author: Titli
"""
import pandas as pd
import glob,os
from itertools import groupby
import multiprocessing 
from joblib import Parallel, delayed
import pickle
import math 
import csv

df= pd.read_excel('Protein ID for TSR Manuscript 8-24-2018.xlsx',sheet_name='PKABC')
group_dict = df.groupby('Name')['PDB IDs'].apply(list).to_dict()
for k,v in group_dict.items():
    print(k, len(v), v)

protList = df['PDB IDs'].tolist()
print(len(protList),protList)

# helper code
def calculate_std(xs):
    mean = sum(xs) / len(xs)   # mean
    var  = sum(pow(x-mean,2) for x in xs) / len(xs)  # variance
    std  = math.sqrt(var)  # standard deviation
    return std

## STEP 1
# find common keys in all proteins
def find_common_keys_allProteins(protList, directory):
    total_keys = []
    total_unq_keys = []
    common_keys = []
 
    for p in protList:
        #print(p)
        df = pd.read_csv(directory+p+".triplets_29_35", sep="\t", usecols=[0],names=['key'])
        total_keys.extend(df['key'].tolist())
        total_unq_keys.extend(list(set(df['key'].tolist())))
    
    total_unq_keys.sort()
    #print(len(total_unq_keys))
    common_unq_keys = [key for key, group in groupby(total_unq_keys) if len(list(group))>=len(protList)]
    #print(len(common_unq_keys))
    
    for p in protList: # get total common keys
        #print(p)
        df = pd.read_csv(directory+p+".triplets_29_35", sep="\t", usecols=[0],names=['key'])
        df_common = df[df['key'].isin(common_unq_keys)]
        common_keys.extend(df_common['key'].tolist())
     
    total_unq_out = open(directory+"allProteins_total_unq_keys.p", "wb")
    common_unq_out = open(directory+"allProteins_common_unq_keys.p", "wb")
    total_out = open(directory+"allProteins_total_keys.p", "wb")
    common_out = open(directory+"allProteins_common_keys.p", "wb")
    pickle.dump(total_unq_keys, total_unq_out)
    pickle.dump(common_unq_keys, common_unq_out)
    pickle.dump(total_keys, total_out)
    pickle.dump(common_keys, common_out)
    
    print("allProteins_total_unq_keys = ",len(total_unq_keys))
    print("allProteins_common_unq_keys = ", len(common_unq_keys))
    print("allProteins_total_keys = ",len(total_keys))
    print("allProteins_common_keys = ", len(common_keys))
    return total_unq_keys, common_unq_keys

# find common keys present in each subgroup which are part of common keys for all proteins
def find_common_keys(grp,protList, directory):
    total_keys = []
    common_keys = []
    common_unq_keys = pickle.load(open(directory+"allProteins_common_unq_keys.p", "rb")) # common unique keys for all proteins
    
    f_out = open(directory+str(grp)+"_commonKey_Details.txt", "w")

    for p in protList:
        #print(p)
        df = pd.read_csv(directory+p+".triplets_29_35", sep="\t", usecols=[0],names=['key'])
        
        df_total = df #with freq
        df_common = df[df['key'].isin(common_unq_keys)]#with freq
        df_uncommon = df[~df['key'].isin(common_unq_keys)]#with freq
    
        f_out.write(str(grp)+"\t"+str(p)+"\t"+str(len(df_total))+"\t"+str(len(df_common))+"\t"+str(len(df_uncommon))+"\n")
        
        total_keys.extend(df_total['key'].tolist())
        common_keys.extend(df_common['key'].tolist())
        
    total_keys.sort()
    common_keys.sort()    #print(len(common_unq_keys))
    uncommon_keys = list(set(total_keys) - set(common_keys))
        
    total_out = open(directory+grp+"_total_keys.p", "wb")
    common_out = open(directory+grp+"_common_keys.p", "wb")
    pickle.dump(total_keys, total_out)
    pickle.dump(common_keys, common_out)
    print("group, total_keys, common_keys, uncommon keys: ")
    print(grp, len(total_keys), len(common_keys), len(uncommon_keys))
    f_out.close()
    return total_keys, common_keys

# find theta and maxdist for total, common, uncommon keys in each subgroup(with frequency)
def find_theta_maxdist(g,pl,directory):
    csvFile = open(directory+str(g)+'_theta_maxdist_details.csv', 'w')
    writer = csv.writer(csvFile)
    
    print(g, len(pl), directory)
    writer.writerow([str(g), str(len(pl)), str(directory)])
    total_keys = pickle.load(open(directory+g+"_total_keys.p", "rb")) # total keys of a subgroup
    common_keys = pickle.load(open(directory+g+"_common_keys.p", "rb")) # common keys for subgroup

    total_unq_keys = list(set(total_keys))
    common_unq_keys = list(set(common_keys))
    #common_unq_keys = pickle.load(open(directory+"allProteins_common_unq_keys.p", "rb")) # common unique keys for all proteins 
    
    print("for this group:total_keys, common_keys: ")
    print(len(total_keys), len(common_keys))
    writer.writerow([str(g), "total_keys", "common_keys"])
    writer.writerow([str(g), str(len(total_keys)), str(len(common_keys))])
    
    print("for this group:total_unq_keys, common_unq_keys: ")
    print(len(total_unq_keys), len(common_unq_keys))
    writer.writerow([str(g), "total_unq_keys", "common_unq_keys"])
    writer.writerow([str(g), str(len(total_unq_keys)), str(len(common_unq_keys))])
    
    theta_total = [] #multiprocessing.Queue()
    theta_common = []
    theta_uncommon = []
    maxDist_total = []
    maxDist_common = []
    maxDist_uncommon = []
    
    no_of_total_keys = 0 #with freq
    no_of_common_keys = 0 #with freq
    no_of_uncommon_keys = 0 #with freq
    
    f_out_k = open(directory+str(g)+"_key_listing.csv", "w")
    f_out_tm = open(directory+str(g)+"_theta_maxdist_listing.csv", "w")
    writer_k = csv.writer(f_out_k)
    writer_tm = csv.writer(f_out_tm)
    
    for p in pl:
        #print(p)
        df = pd.read_csv(directory+p+".triplets_29_35", sep="\t", usecols=[0,8,10],names=['key','theta','maxDist'])
        df_c = df.loc[df['key'].isin(common_unq_keys)]
        df_u = df.loc[~df['key'].isin(common_unq_keys)]
        df_c.to_csv(directory+p+".common_triplets_29_35", sep="\t", index=False)
        df_u.to_csv(directory+p+".uncommon_triplets_29_35", sep="\t", index=False)

        no_of_total_keys += len(df)
        no_of_common_keys += len(df_c)
        no_of_uncommon_keys += len(df_u)
        writer_k.writerow([str(g), str(p), str(no_of_total_keys), str(no_of_common_keys), str(no_of_uncommon_keys)])
        #print(p, "   no_of_total_keys, no_of_common_keys, no_of_uncommon_keys\n",no_of_total_keys, no_of_common_keys, no_of_uncommon_keys) 

        theta_total.extend(df['theta'].tolist())
        theta_common.extend(df_c['theta'].tolist())
        theta_uncommon.extend(df_u['theta'].tolist())
        maxDist_total.extend(df['maxDist'].tolist())
        maxDist_common.extend(df_c['maxDist'].tolist())
        maxDist_uncommon.extend(df_u['maxDist'].tolist())
        writer_tm.writerow([str(g), str(p), str(len(theta_total)), str(len(theta_common)), str(len(theta_uncommon)), str(len(maxDist_total)), str(len(maxDist_common)), str(len(maxDist_uncommon))])

    #print(len(theta_total),len(theta_common),len(theta_uncommon),len(maxDist_total),len(maxDist_common),len(maxDist_uncommon))
    avg_theta_total = round(sum(theta_total)/len(theta_total),2)
    avg_theta_common = round(sum(theta_common)/len(theta_common),2)
    avg_theta_uncommon = round(sum(theta_uncommon)/len(theta_uncommon),2)
    avg_maxDist_total = round(sum(maxDist_total)/len(maxDist_total),2)
    avg_maxDist_common = round(sum(maxDist_common)/len(maxDist_common),2)
    avg_maxDist_uncommon = round(sum(maxDist_uncommon)/len(maxDist_uncommon),2)

    std_theta_total = round(calculate_std(theta_total),4)
    std_theta_common = round(calculate_std(theta_common),4)
    std_theta_uncommon = round(calculate_std(theta_uncommon),4)

    std_maxDist_total = round(calculate_std(maxDist_total),4)
    std_maxDist_common = round(calculate_std(maxDist_common),4)
    std_maxDist_uncommon = round(calculate_std(maxDist_uncommon),4)
    
    print("Theta: ", avg_theta_total, avg_theta_common, avg_theta_uncommon, avg_maxDist_total, avg_maxDist_common, avg_maxDist_uncommon)
    print("std ", std_theta_total, std_theta_common, std_theta_uncommon, std_maxDist_total, std_maxDist_common, std_maxDist_uncommon)
    print("no_of_total_keys, no_of_common_keys, no_of_uncommon_keys:")
    print(no_of_total_keys, no_of_common_keys, no_of_uncommon_keys)
    writer.writerow(["Theta: ", avg_theta_total, avg_theta_common, avg_theta_uncommon, avg_maxDist_total, avg_maxDist_common, avg_maxDist_uncommon])
    writer.writerow(["std ", std_theta_total, std_theta_common, std_theta_uncommon, std_maxDist_total, std_maxDist_common, std_maxDist_uncommon])
    writer.writerow(["no_of_total_keys", "no_of_common_keys", "no_of_uncommon_keys:"])
    writer.writerow([no_of_total_keys, no_of_common_keys, no_of_uncommon_keys])
    
    f_out_k.close()
    f_out_tm.close()
    csvFile.close()
    #return(avg_theta_total, avg_theta_common, avg_theta_uncommon, avg_maxDist_total, avg_maxDist_common, avg_maxDist_uncommon)    
     
num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)

directory = "./../Dataset/lexiographic/"

# func call: step 1
print("\n\nFinding common keys from all proteins..........")
find_common_keys_allProteins(protList, directory)

#func call: step 2
print ("\n\nFinding common keys for each subgroup..........") 
#Parallel(n_jobs=num_cores, verbose=50)(delayed(find_common_keys(g,pl, directory))(pl)for g,pl in group_dict.items())
for  g,pl in group_dict.items():
    find_common_keys(g,pl, directory)

# func call: step 3
print("\n\nFinding theta maxdist ............")
#Parallel(n_jobs=num_cores, verbose=50)(delayed(find_theta_maxdist(g,pl, directory))(pl)for g,pl in group_dict.items())
for g,pl in group_dict.items():
    find_theta_maxdist(g,pl, "./../Dataset/lexiographic/")
