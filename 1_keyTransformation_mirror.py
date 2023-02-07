# -*- coding: utf-8 -*-
'''
Code for calculating PDB to TSRs with changed length bins (domain knowledge incorporated)
'''
import math
import multiprocessing
from joblib import Parallel, delayed
import pandas as pd
import os, argparse, time
import numpy as np

start_time = time.time()

dTheta = 30
dLen = 35
numOfLabels = 20  # Without Amino Acid grouping

data_dir="./../Dataset/" # output dir
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
subdir="lexiographic/mirror/"
if not os.path.exists(data_dir+subdir):
    os.makedirs(data_dir+subdir)

#files= [f.split("/")[3].split(".")[0] for f in glob.glob(data_dir+"*.pdb")]
#print(files)
#print(len(files))

parser = argparse.ArgumentParser(description='Key Calculation with chain info + mirror.')
parser.add_argument('sample_details', default='sample_details.csv',type=str, help='Enter input file for protein list to be downloaded. File format should be .csv file with first column having list of proteins with header=protein')
args = parser.parse_args()

df= pd.read_csv(args.sample_details,sep=',',header=0)
df['protchain'] = df['protein']+df['chain']
files = df['protchain'].tolist()
print(len(files))
#chain_dict = df.groupby('chain')['protein'].apply(list).to_dict()
#for k,v in chain_dict.items():
    #print (k,v)

def thetaClass_(Theta):
    #classT=0
    if Theta<=0.0001:
        classT=30
    elif Theta>0.0001 and Theta<12.11:
        classT=1
    elif Theta>=12.11 and Theta<17.32:
        classT=2
    elif Theta>=17.32 and Theta<21.53:
        classT=3
    elif Theta>=21.53 and Theta<25.21:
        classT=4
    elif Theta>=25.21 and Theta<28.54:
        classT=5
    elif Theta>=28.54 and Theta<31.64:
        classT=6
    elif Theta>=31.64 and Theta<34.55:
        classT=7
    elif Theta>=34.55 and Theta<37.34:
        classT=8
    elif Theta>=37.34 and Theta<40.03:
        classT=9
    elif Theta>=40.03 and Theta<42.64:
        classT=10
    elif Theta>=42.64 and Theta<45.17:
        classT=11
    elif Theta>=45.17 and Theta<47.64:
        classT=12
    elif Theta>=47.64 and Theta<50.05:
        classT=13
    elif Theta>=50.05 and Theta<52.43:
        classT=14
    elif Theta>=52.43 and Theta<54.77:
        classT=15
    elif Theta>=54.77 and Theta<57.08:
        classT=16
    elif Theta>=57.08 and Theta<59.38:
        classT=17
    elif Theta>=59.38 and Theta<61.64:
        classT=18
    elif Theta>=61.64 and Theta<63.87:
        classT=19
    elif Theta>=63.87 and Theta<66.09:
        classT=20
    elif Theta>=66.09 and Theta<68.30:
        classT=21
    elif Theta>=68.30 and Theta<70.5:
        classT=22
    elif Theta>=70.5 and Theta<72.69:
        classT=23
    elif Theta>=72.69 and Theta<79.2:
        classT=24
    elif Theta>=79.2 and Theta<81.36:
        classT=25
    elif Theta>=81.36 and Theta<83.51:
        classT=26
    elif Theta>=83.51 and Theta<85.67:
        classT=27
    elif Theta>=85.67 and Theta<87.80:
        classT=28
    elif Theta>=87.80 and Theta<=90.00:
        classT=29
    return classT

def dist12Class_(dist12):
    #classL=0
    if (dist12<3.83):
        classL=1
    elif dist12>=3.83 and dist12<7.00:
        classL=2
    elif dist12>=7.00 and dist12<9.00:
        classL=3
    elif dist12>=9.00 and dist12<11.00:
        classL=4
    elif dist12>=11.00 and dist12<14.00:
        classL=5
    elif dist12>=14.00 and dist12<17.99:
        classL=6
    elif dist12>=17.99 and dist12<21.25:
        classL=7
    elif dist12>=21.25 and dist12<23.19:
        classL=8
    elif dist12>=23.19 and dist12<24.8:
        classL=9
    elif dist12>=24.8 and dist12<26.26:
        classL=10
    elif dist12>=26.26 and dist12<27.72:
        classL=11
    elif dist12>=27.72 and dist12<28.9:
        classL=12
    elif dist12>=28.9 and dist12<30.36:
        classL=13
    elif dist12>=30.36 and dist12<31.62:
        classL=14
    elif dist12>=31.62 and dist12<32.76:
        classL=15
    elif dist12>=32.76 and dist12<33.84:
        classL=16
    elif dist12>=33.84 and dist12<35.13:
        classL=17
    elif dist12>=35.13 and dist12<36.26:
        classL=18
    elif dist12>=36.26 and dist12<37.62:
        classL=19
    elif dist12>=37.62 and dist12<38.73:
        classL=20
    elif dist12>=38.73 and dist12<40.12:
        classL=21
    elif dist12>=40.12 and dist12<41.8:
        classL=22
    elif dist12>=41.8 and dist12<43.41:
        classL=23
    elif dist12>=43.41 and dist12<45.55:
        classL=24
    elif dist12>=45.55 and dist12<47.46:
        classL=25
    elif dist12>=47.46 and dist12<49.69:
        classL=26
    elif dist12>=49.69 and dist12<52.65:
        classL=27
    elif dist12>=52.65 and dist12<55.81:
        classL=28
    elif dist12>=55.81 and dist12<60.2:
        classL=29
    elif dist12>=60.2 and dist12<64.63:
        classL=30
    elif dist12>=64.63 and dist12<70.04:
        classL=31
    elif dist12>=70.04 and dist12<76.15:
        classL=32
    elif dist12>=76.15 and dist12<83.26:
        classL=33
    elif dist12>=83.26 and dist12<132.45:
        classL=34
    elif dist12>=132.45:
        classL=35
    return classL

def calcDist(indexLabel1,indexLabel2):
    x1=xCord[indexLabel1]
    x2=xCord[indexLabel2]
    y1=yCord[indexLabel1]
    y2=yCord[indexLabel2]
    z1=zCord[indexLabel1]
    z2=zCord[indexLabel2]
    distance=(((x1-x2)**2+(y2-y1)**2+(z2-z1)**2)**0.5)
    return distance

def indexFind(index_of_2,i1,j1,k1):
    if index_of_2==i1:
        indexOf0=j1
        indexOf1=k1
    elif index_of_2==j1:
        indexOf0=i1
        indexOf1=k1
    elif index_of_2==k1:
        indexOf0=i1
        indexOf1=j1

    return indexOf0, indexOf1

aminoAcidCode=open("aminoAcidCode_lexicographic_new.txt","r")

aminoAcidLabel={}
#aminoAcidGroup={}
for amino in aminoAcidCode:
    amino=amino.split()
    aminoAcidLabel[amino[0]]=int(amino[1])
aminoAcidCode.close()

#for fileName in files:
def parallelcode(fileName): #fileName includes chain name
    protName = fileName[:-1]
    chain = fileName[-1].upper()
    print (fileName, protName, chain)
    #start_time=time.time()
    filesDict={}
    inFile=open(data_dir+protName+'.pdb','r')
    outFile2 = open(data_dir+subdir+fileName+".keys_theta30_maxdist35_mirror", "w")
    fileTriplets = open(data_dir+subdir+fileName+".triplets_theta30_maxdist35_mirror","w")
    # write header 
    fileTriplets.writelines("key\taa0\tpos0\taa1\tpos1\taa2\tpos2\tclassT1\tTheta\tclassL1\tmaxDist\tx0\ty0\tz0\tx1\ty1\tz1\tx2\ty2\tz2\tTheta1\n")
    
    global xCord, yCord, zCord
    aminoAcidName={}
    xCord={}
    yCord={}
    zCord={}
    seq_number={}
    counter=0
    flag=False
    for i in inFile:
        if flag==True:
            break
        #if (i[0:6].rstrip()=="NUMMDL"):
            #numOfModels=i[10:14].rstrip()
        if ((i[0:6].rstrip()=="ENDMDL") or (i[0:6].rstrip()=='TER' and i[21].rstrip()==chain)):
            flag=True
            break
        if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
            break
         
        if(i[0:4].rstrip())=="ATOM"and(i[13:15].rstrip())=="CA"and(i[16]=='A' or i[16]==' ') and i[21:22].strip()==chain and i[17:20]!= "UNK" :
            #print (i)
            aminoAcidName[counter]=int(aminoAcidLabel[i[17:20]])
            xCord[counter]=(float(i[30:38]))
            yCord[counter]=(float(i[38:46]))
            zCord[counter]=(float(i[46:54]))
            seq_number[counter]=str(i[22:27])
            counter+=1
            
    protLen=len(yCord)
    initialLabel=[]
    sortedLabel=[]
    sortedIndex=[]
    for m in range(0,3):
        initialLabel.append(0)
        sortedLabel.append(0)
        sortedIndex.append(0)
    for i in range(0,protLen-2):
        for j in range(i+1,protLen-1):
            for k in range(j+1, protLen):
                global i1,j1,k1
                i1=i
                j1=j
                k1=k
                keepLabelIndex={}
                keepLabelIndex[aminoAcidName[i]]=i
                keepLabelIndex[aminoAcidName[j]]=j
                keepLabelIndex[aminoAcidName[k]]=k
                initialLabel[0]=aminoAcidName[i]
                initialLabel[1]=aminoAcidName[j]
                initialLabel[2]=aminoAcidName[k]
                sortedLabel=list(initialLabel)
                sortedLabel.sort(reverse=True)
                if (sortedLabel[0]==sortedLabel[1])and(sortedLabel[1]==sortedLabel[2]):
                    dist1_2Temp=calcDist(i,j)
                    dist1_3Temp=calcDist(i,k)
                    dist2_3Temp=calcDist(j,k)
                    if dist1_2Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                        indexOf0=i
                        indexOf1=j
                        indexOf2=k
                    elif dist1_3Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                        indexOf0=i
                        indexOf1=k
                        indexOf2=j
                    else:
                        indexOf0=j
                        indexOf1=k
                        indexOf2=i
                elif(aminoAcidName[i]!=aminoAcidName[j])and(aminoAcidName[i]!=aminoAcidName[k])and(aminoAcidName[j]!=aminoAcidName[k]):
                    for index_ in range(0,3):
                        sortedIndex[index_]=keepLabelIndex[sortedLabel[index_]]
                    indexOf0=sortedIndex[0]
                    indexOf1=sortedIndex[1]
                    indexOf2=sortedIndex[2]
                elif(sortedLabel[0]==sortedLabel[1])and(sortedLabel[1]!=sortedLabel[2]):
                    indexOf2=keepLabelIndex[sortedLabel[2]]
                    indices=indexFind(indexOf2,i,j,k)
                    a=indexOf2
                    b=indices[0]
                    c=indices[1]
                    dist1_3Temp=calcDist(b,a)
                    dist2_3Temp=calcDist(c,a)
                    if dist1_3Temp>=dist2_3Temp:
                        indexOf0=indices[0]
                        indexOf1=indices[1]	
                    else:
                        indexOf0=indices[1]
                        indexOf1=indices[0]
                elif(sortedLabel[0]!=sortedLabel[1])and(sortedLabel[1]==sortedLabel[2]):
                    indexOf0=keepLabelIndex[sortedLabel[0]]
                    indices=indexFind(indexOf0,i,j,k)
                    if calcDist(indexOf0,indices[0])>= calcDist(indexOf0,indices[1]):
                        indexOf1=indices[0]
                        indexOf2=indices[1]	
                    else:
                        indexOf2=indices[0]
                        indexOf1=indices[1]
                dist01=calcDist(indexOf0,indexOf1)
                s2=dist01/2
                dist02=calcDist(indexOf0,indexOf2)
                s1=dist02
                #dist12=dist01
                dist03=calcDist(indexOf1,indexOf2)
                maxDist=max(dist01,dist02,dist03)
                s3=(((xCord[indexOf0]+xCord[indexOf1])/2-xCord[indexOf2])**2+((yCord[indexOf0]+yCord[indexOf1])/2-yCord[indexOf2])**2+((zCord[indexOf0]+zCord[indexOf1])/2-zCord[indexOf2])**2)**0.5
                Theta1=180*(math.acos((s1**2-s2**2-s3**2)/(2*s2*s3)))/3.14
                if Theta1<=90:
                    Theta=Theta1
                else:
                    Theta=abs(180-Theta1)
                classT1=thetaClass_(Theta)
                classL1=dist12Class_(maxDist)

                ##getting the positions of AminoAcids in sequence
                position0 = str(list(seq_number.values())[indexOf0])
                position1 = str(list(seq_number.values())[indexOf1])
                position2 = str(list(seq_number.values())[indexOf2])

                aacd0 = list(aminoAcidLabel.keys())[list(aminoAcidLabel.values()).index(aminoAcidName[indexOf0])]
                aacd1 = list(aminoAcidLabel.keys())[list(aminoAcidLabel.values()).index(aminoAcidName[indexOf1])]
                aacd2 = list(aminoAcidLabel.keys())[list(aminoAcidLabel.values()).index(aminoAcidName[indexOf2])]

                x0 = str(xCord.get(indexOf0))
                y0 = str(yCord.get(indexOf0))
                z0 = str(zCord.get(indexOf0))

                x1 = str(xCord.get(indexOf1))
                y1 = str(yCord.get(indexOf1))
                z1 = str(zCord.get(indexOf1))

                x2 = str(xCord.get(indexOf2))
                y2 = str(yCord.get(indexOf2))
                z2 = str(zCord.get(indexOf2))

                key_2=dLen*dTheta*(numOfLabels**2)*(aminoAcidName[indexOf0]-1)+\
                      dLen*dTheta*(numOfLabels)*(aminoAcidName[indexOf1]-1)+\
                      dLen*dTheta*(aminoAcidName[indexOf2]-1)+\
                          dTheta*(classL1-1)+(classT1-1)
                
                # negetion of keys for mirror      
                if (Theta1 > 90) & (aminoAcidName[indexOf0] > aminoAcidName[indexOf1]):  #& (aminoAcidName[indexOf0] != aminoAcidName[indexOf1])
                    key_2 = (-1)*key_2
                #print (key_2)
                                
                if key_2 in filesDict:
                    filesDict[key_2]+=1
                else:
                    filesDict[key_2]=1

                line = (str(key_2)+"\t"+\
                        str(aacd0)+"\t"+str(position0)+"\t"+str(aacd1)+"\t"+str(position1)+"\t"+str(aacd2)+"\t"+str(position2)+"\t"+\
                            str(classT1)+"\t"+str(Theta)+"\t"+str(classL1)+"\t"+str(maxDist)+"\t"+\
                            x0+"\t"+y0+"\t"+z0+"\t"+x1+"\t"+y1+"\t"+z1+"\t"+x2+"\t"+y2+"\t"+z2+"\t"+\
                            str(Theta1)+"\n")
                fileTriplets.writelines(line)
                #print (line)

    for value_ in filesDict:
        #print (value_)
        outFile2.writelines([str(value_),'\t', str(filesDict[value_]),'\n'])

    print ("FILENAME=",fileName,'\t',"NUM OF AMINOACIDS=",protLen)
    outFile2.close()
    fileTriplets.close()
    ## end of parallelcode()

num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)
#for chainId,files in chain_dict.items():
#Parallel(n_jobs=num_cores, verbose=50)(delayed(parallelcode)(fileName,chainId.upper())for fileName in files)
Parallel(n_jobs=num_cores, verbose=50)(delayed(parallelcode)(fileName)for fileName in files)

print("CODE END\t TIME=(min)", np.round((time.time()-start_time)/60,2))

