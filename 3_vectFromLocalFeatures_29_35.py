import copy,glob,os
import pandas as pd

subdir = './../Dataset/lexiographic/'
#filesList= [f.split("/")[3].split(".")[0] for f in glob.glob("./../Dataset/*.pdb")]
df= pd.read_excel('Protein ID for TSR Manuscript 8-24-2018.xlsx',sheet_name='PKABC')
filesList = df['PDB IDs'].tolist()
print(len(filesList), filesList)

#saves all the counts for each key 
allKeyCountDict={}

f1_in=open(subdir+'localFeatureSelection_44PKABC_29_35.txt','r')
f1_out=open(subdir+'localFeatureVect_44PKABC_29_35.csv','w')

keyDict={}
for i in f1_in:
    i=i.split()
    keyDict[i[0].rstrip()]=0
print (len(keyDict))

for i in filesList:
    f2_in=open(subdir+i+'.keys_29_35','r')
    keyDict1=copy.deepcopy(keyDict)
    
    for j in f2_in:
    	j=j.split()
    	if j[0].rstrip() in keyDict1:
    	    keyDict1[j[0].rstrip()]=j[1].rstrip()
    for k in keyDict1:
        f1_out.writelines([str(keyDict1[k]),','])
    f1_out.writelines(['0','\n'])
    	
print ("Code End.")
