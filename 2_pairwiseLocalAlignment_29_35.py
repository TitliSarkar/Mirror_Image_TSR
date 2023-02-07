import math,glob
import numpy
import os 
import pandas as pd

subdir = './../Dataset/lexiographic/'
#filesList= [f.split("/")[3].split(".")[0] for f in glob.glob("./../Dataset/*.pdb")]
df= pd.read_excel('Protein ID for TSR Manuscript 8-24-2018.xlsx',sheet_name='PKABC')
filesList = df['PDB IDs'].tolist()
print(len(filesList), filesList)

#saves all the counts for each key 
allKeyCountDict={}
n = len(filesList)

for i in filesList:
    f1=open(subdir+i+'.keys_29_35','r')
    for j in f1:
        j=j.split()
        key_=j[0].rstrip()
        val_=int(j[1].rstrip())
        if key_ in allKeyCountDict:
            allKeyCountDict[key_].append(val_)
        else:
            allKeyCountDict[key_]=[]
            allKeyCountDict[key_].append(val_)

f2_out=open(subdir+'localFeatureSelection_44PKABC_29_35.txt','w')

for keys_ in allKeyCountDict:
    list_=[]
    list_=allKeyCountDict[keys_]
    numOfMatch=len(list_) #number of proteins that have the key
    numOfGap=n-numOfMatch #number of proteins that DO NOT have the key
    mean_=numpy.mean(list_) #average number of times a key occurs in the list of proteins
    mad_=0.0
   
    for i in range(len(list_)):
        mad_=mad_+math.fabs(list_[i]-mean_)
    mad_=float(mad_)/float(numOfMatch)# calculate MAD
 
    #if ((numOfMatch<=math.ceil(n/5) and mad_<=0.01)):# and (numOfMatch>2)):
    f2_out.writelines([str(keys_),'\t',str(mad_),'\n'])
  
print ("Code End.")
