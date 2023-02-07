import os, glob
#import itertools
import time
import math
import numpy as np 
import pandas as pd
import multiprocessing 
from joblib import Parallel, delayed

dataset = os.getcwd().split("/")[-2]
print(dataset)

directory = "./../Dataset/lexiographic/"

df= pd.read_excel('Protein ID for TSR Manuscript 8-24-2018.xlsx',sheetname='PKABC')
chain_dict =  dict(zip(df['PDB IDs'], df['Chain']))
for k,v in chain_dict.items():
    print (k,v)

def countAminos(p):
    chain = chain_dict[p].upper()
    inFile=open("./../Dataset/"+p+'.pdb','r')
    counter = 0

    for i in inFile:
        if (i[0:6].rstrip()=="NUMMDL"):
            numOfModels=i[10:14].rstrip()
        if ((i[0:6].rstrip()=="ENDMDL") or (i[0:6].rstrip()=='TER' and i[21].rstrip()==chain)):
            break
        if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
            break
         
        if(i[0:4].rstrip())=="ATOM"and(i[13:15].rstrip())=="CA"and(i[16]=='A'or i[16]==' ') and i[21:22].strip()==chain and i[17:20]!= "UNK" :
            counter+=1
    return counter

# get amino acid codes in a dict: label:[list of Amino Names]
amino_label_dict = {}
with open("aminoAcidCode_grouping_type0.txt", "r") as f:
    for line in f:
        if not line.strip():
            continue
        line = line.strip().split()
        if line[1] in amino_label_dict:
            amino_label_dict[line[1]].append(line[0])
        else:
            amino_label_dict[line[1]]=[]
            amino_label_dict[line[1]].append(line[0])
f.close()

#seperate unchanged and changed aminos
unchanged = [] # list of unchanged aminos
changed = [] # list of all changed aminos 
for k,v in amino_label_dict.items():
    #print(k,v)
    if len(v)==1:
        unchanged.append(v)
    else:
        changed.append(v)
unchanged = [y for x in unchanged for y in x] 
unchanged.sort() 
changed = [y for x in changed for y in x] 
changed.sort() 
print("Singleton aminos are : ",unchanged) #['CYS', 'GLY', 'HIS', 'MET', 'PRO', 'TYR']
print("Grouped aminos are : ", changed) # [['ASN', 'GLN'], ['SER', 'THR'], ['LYS', 'ARG'], ['ALA', 'VAL'], ['LEU', 'ILE'], ['PHE', 'TRP'], ['ASP', 'GLU']]

# get the list of the proteins
filesList = [f.split("/")[4].split(".")[0] for f in glob.glob(directory+"*.triplets_29_35")]
print("#of lex files: ",len(filesList), filesList)

def calculate_TSR_percentage_parallel(p, unchanged, changed):
    print(p)
    count_total = 0
    count_u_u_u = 0
    count_u_u_c = 0
    count_u_c_c = 0
    count_c_c_c = 0
    
    f_in= open(directory+p+".triplets_29_35","r")
    f_out= open(directory+p+"_TSR_percentage_temp.txt","w")

    for lines in f_in:
        if not lines.strip():
            continue
        line = lines.strip().split()
        if len(line) != 20:
            continue
        count_total += 1
        
        triplets = [line[1], line[3], line[5]]
        status = []
        for amino in triplets:
            if amino in unchanged:
                status.append('u')
            if amino in changed:
                status.append('c')
        #print(status)
        if status.count('u') == 3:
            count_u_u_u += 1
            print(status, count_u_u_u)
        if status.count('u') == 2 and status.count('c') == 1:
            count_u_u_c += 1
            print(status, count_u_u_c)
        if status.count('u') == 1 and status.count('c') == 2:
            count_u_c_c += 1         
            print(status, count_u_c_c)
        if status.count('c') == 3:
            count_c_c_c += 1         
            print(status, count_c_c_c)

    
    #total_singleton_TSR = count_u_u_u
    #total_grouped_TSR  = (count_u_u_c + count_u_c_c + count_c_c_c)
    
    #total_singleton_TSR_percentage = np.round((total_singleton_TSR / count_total), 2)
    #total_grouped_TSR_percentage = np.round((total_grouped_TSR / count_total), 2)
    
    #ratio_total = np.round((total_grouped_TSR/total_singleton_TSR),2)
    
    line = [str(p)]
    line.append(str(countAminos(p)))
    line.append(str(count_total))
    line.append(str(count_u_u_u))
    line.append(str(count_u_u_c))
    line.append(str(count_u_c_c))
    line.append(str(count_c_c_c))
    line.append("\n")
    print(line)
    f_out.write("\t".join(line))
    f_out.close()
    f_in.close()
    #end of func

# get the number of cores in multiprocessing
num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)

Parallel(n_jobs=num_cores, verbose=50)(delayed(calculate_TSR_percentage_parallel)(prot, unchanged, changed)for prot in filesList)
# concatinate output files
f_out = open(directory+dataset+"_TSR_percentage.txt", "w")
header = ['protein','#of amino acids','total','total u_u_u','total u_u_c','total u_c_c','total c_c_c','\n']
f_out.write("\t".join(header))
for f in filesList:
    f_in = open(directory+f+"_TSR_percentage_temp.txt", "r")
    for line in f_in:
        f_out.write(line)
    f_in.close()
    os.remove(directory+f+"_TSR_percentage_temp.txt")
f_out.close()
print("Code End.")
