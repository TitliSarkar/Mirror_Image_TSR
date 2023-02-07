import os, glob
#import itertools
import time
import math
import pandas as pd
import multiprocessing 
from joblib import Parallel, delayed

directory = "./../Dataset/lexiographic/"
os.makedirs(directory, exist_ok=True)

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
print("Unchanged aminos are : ",unchanged) #['CYS', 'GLY', 'HIS', 'MET', 'PRO', 'TYR']
print("Changed aminos are : ", changed) # [['ASN', 'GLN'], ['SER', 'THR'], ['LYS', 'ARG'], ['ALA', 'VAL'], ['LEU', 'ILE'], ['PHE', 'TRP'], ['ASP', 'GLU']]

# get the list of the proteins
filesList = [f.split("/")[4].split(".")[0] for f in glob.glob(directory+"*.triplets_29_35")]
print("#of lex files: ",len(filesList), filesList)

def calculate_amino_acid_percentage_parallel(p, unchanged, changed):
    print(p)
    count_total = 0
    count_changed_dict = {}
    count_unchanged_dict = {}
    for c in changed:
        count_changed_dict[c] = 0
    for u in unchanged:
        count_unchanged_dict[u] = 0

    f_in= open(directory+p+".triplets_29_35","r")
    f_out= open(directory+p+"_amino_percentage_temp.txt","w")
    
    for lines in f_in:
        if not lines.strip():
            continue
        line = lines.strip().split()
        
        count_total += 3
        
        if line[1] in unchanged:
            count_unchanged_dict[line[1]] += 1
        if line[3] in unchanged:
            count_unchanged_dict[line[3]] += 1
        if line[5] in unchanged:
            count_unchanged_dict[line[5]] += 1 
     
        if line[1] in changed:
            count_changed_dict[line[1]] += 1
        if line[3] in changed:
            count_changed_dict[line[3]] += 1
        if line[5] in changed:
            count_changed_dict[line[5]] += 1
    
    #print(count_unchanged_dict)
    #print(count_changed_dict)
    tmp_dict = {}
    total_unchanged_percentage = 0.00
    total_changed_percentage = 0.00

    for key in sorted(count_unchanged_dict):
        tmp_dict[key] = count_unchanged_dict[key]
        per_ = round((count_unchanged_dict[key] / count_total)*100, 2)
        total_unchanged_percentage += per_
    for key in sorted(count_changed_dict):
        tmp_dict[key] = count_changed_dict[key]
        per_ = round((count_changed_dict[key] / count_total)*100, 2)
        total_changed_percentage += per_


    line = [str(p)]
    line.append("#of amino acids")
    line.append(str(countAminos(p)))
    line.append("total")
    line.append(str(count_total))
    line.append("total_unchaged_percetage")
    line.append(str(total_unchanged_percentage))
    line.append("total_chaged_percetage")
    line.append(str(total_changed_percentage))
    for key in sorted(tmp_dict):
        line.append(str(key))
        line.append(str(tmp_dict[key]))
        #line.append(str(round((tmp_dict[key] / count_total)*100), 2))
        per_ = round((tmp_dict[key] / count_total)*100, 2)
        line.append(str(per_))
    line.append("\n")
   
    print(line)
    f_out.write("\t".join(line))
    f_out.close()
    f_in.close()
    #end of func

# get the number of cores in multiprocessing
num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)

Parallel(n_jobs=num_cores, verbose=50)(delayed(calculate_amino_acid_percentage_parallel)(prot, unchanged, changed)for prot in filesList)

# concatinate output files
f_out = open(directory+"44PKABC_amino_percentage.txt", "w")
for f in filesList:
    f_in = open(directory+f+"_amino_percentage_temp.txt", "r")
    for line in f_in:
        f_out.write(line)
    f_in.close()
    os.remove(directory+f+"_amino_percentage_temp.txt")
f_out.close()
print("Code End.")
