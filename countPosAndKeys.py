#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 11:45:32 2018

@author: txs7980
"""
import os,glob
import pandas as pd

subdir = './../Dataset/lexiographic/'

filesList= sorted([y.split("/")[3].split(".")[0] for y in glob.glob("./../Dataset/*.pdb")])
print(len(filesList), filesList)

def countAminos(chain_dict):
    amino_dict ={}
    for chainId,files in chain_dict.items():
        for fileName in files:
            inFile=open(fileName+'.pdb','r')
            counter=0
            for i in inFile:
                #if (i[0:6].rstrip()=="NUMMDL"):
                    #numOfModels=i[10:14].rstrip()
                if ((i[0:6].rstrip()=="ENDMDL") or (i[0:6].rstrip()=='TER' and i[21].rstrip()==chainId)):
                    break
                
                if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
                    break
                
                if(i[0:4].rstrip())=="ATOM"and(i[13:15].rstrip())=="CA"and(i[16]=='A'or i[16]==' ' and i[21:22]==chainId)and i[17:20]!= "UNK" :
                    counter+=1
            amino_dict[fileName] = int(counter)
    for i in filesList:
        print(i, amino_dict[i])

def countKeys():
    for f in filesList:
        inFile=open(f+'.commonTriplets_29_35_maxdist11','r') #commonTriplets_PKB_notOthers_29_35
        keys=0
        for i in inFile:
            #print(i)
            #counter +=int(i.split()[1]) #total keys #i.split[1] in case of .keys 
            keys +=1 
        print (keys)
        
def count_sum(): # count total #connected keys from scalaprocessed file
    subdir = '37RE_21_12_freq1\\' 
    os.chdir(subdir)
    print(os.getcwd())
    
    for i in filesList:
        df = pd.read_csv(i+'AllCCs.csv')
        print(df['size_of_component'].sum())

'''def search_a_key_in_every_protein(fileList, key): # search a single key
    present = []
    for f in fileList:
        df = pd.read_csv(f+'.triplets_29_35',sep="\t",names=['key','freq'])
       # df = pd.read_csv(f+'.triplets_29_35', sep='\t', names = ['key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'])

        keylist = list(set(df['key'].tolist()))
        if key in keylist:
            print (f, "present")
            present.append(f)
        else:
            print(f)
    print(len(present),len(fileList))'''
        
def search_a_key_in_every_protein_return_valid_keys(fileList, keys): # searching if all/some of the keys are present each protein
    print("...............................................>")
    present = []
    present_keys = []
    #present_df = pd.DataFrame(columns=['key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'])
    #i = 1
    for f in fileList:
        print(f)
        df = pd.read_csv(subdir+f+'.triplets_29_35', sep='\t', names = ['key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'])
        keylist = list(set(df['key'].tolist()))
        if len(set(keys) & set(keylist)) != 0:
            present.append(f)
            #print(df.loc[df['key'] == list(set(keys) & set(keylist))[0]])
            #present_df.append(df.loc[df['key'] == list(set(keys) & set(keylist))[0]])
            #i += 1
            present_keys.extend(list(set(keys) & set(keylist)))
    if (len(present) == len(fileList)):
        print("present", present_keys)
        #print(present_df)
        return True,present_keys
    else:
        print(" not present !", present_keys)
        return False,_

def search_a_key_in_every_protein(fileList, keys): # searching if all/some of the keys are present each protein
    present = []
    for f in fileList:
        print(f)
        df = pd.read_csv(f+'.keys_29_35',sep="\t",names=['key','freq'])
        keylist = list(set(df['key'].tolist()))
        if len(set(keys) & set(keylist)) != 0:
            present.append(f)
            
    presentdict = {} # among the present proteins, how many proteins from different groups
    for p in present:
        if kinaseDict[p] in presentdict:
            presentdict[kinaseDict[p]].append(p)
        else:
            presentdict[kinaseDict[p]] = []
            presentdict[kinaseDict[p]].append(p)
    for k,v in presentdict.items():
        print (k, len(v), '\n')
    
    print("present: ", len(present), "fileList: ", len(fileList))

def search_combined_key_in_every_protein(fileList, keyList): #search if two keys are present together, each with own stretched range
    for i in range(3):
        present = [] #list of proteins where criteria is satisfied
        
        for f in fileList: # for each prot, search each key in its range
            #print(f)
            df = pd.read_csv(subdir+f+'.keys_29_35',sep="\t",names=['key','freq'])
            keylist = list(set(df['key'].tolist()))
            
            count = 0 # holds number of searched keys present
            for key in keyList:
                stretched_key = list(range(int(key)-i, int(key)+i+1))
                #stretched_key = [int(key)]
                if len(set(keylist) & set(stretched_key)) != 0:
                    count += 1
                '''if key == 5406015:
                    stretched_key = list(range(int(key)-2, int(key)+2+1))
                    if len(set(keylist) & set(stretched_key)) != 0:
                        count += 1
                if key == 5548103:
                    stretched_key = list(range(int(key)-5, int(key)+5+1))
                    if len(set(keylist) & set(stretched_key)) != 0:
                        count += 1'''
            #print(count)
            if count == len(keyList):
                #print (f, '--->  ', kinaseDict[str(f)])
                present.append(f)
        '''presentdict = {} # among the present proteins, how many proteins from different groups
        for p in present:
            if kinaseDict[p] in presentdict:
                presentdict[kinaseDict[p]].append(p)
            else:
                presentdict[kinaseDict[p]] = []
                presentdict[kinaseDict[p]].append(p)
        #for k,v in presentdict.items():
            #print (k, len(v), '\n')'''
		
        print("fileList, present, +/= ")
        print(len(fileList), len(present), i)
    
##### function call 
#countAminos(chain_dict)
#countKeys()
#count_sum()
#key = '2988275' # 4273260 4679259 4679260 5511589             9209176 8803162
#stretched_key = list(range(int(key)-1, int(key)+1+1)) # do +/-10 for each key
#search_a_key_in_every_protein_return_valid_keys(fileList, [])
#valid_stretched_key = [7944448, 7944449, 7944450, 7944451, 7944452, 7944453, 7944454, 7944455, 7944456, 7944457, 7944458, 7944459, 7944460, 7944461, 7944462]
#search_a_key_in_every_protein(filesList,stretched_key) #['7903915', '3803315']
#print(subdir, key)
#print ("Code End")
#print(subdir, key)

search_combined_key_in_every_protein(filesList,[6362130]) 
