# -*- coding: utf-8 -*-
''' In Use '''
import os,glob
import pandas as pd

subdirName = './../Dataset/lexiographic/'

filesList= sorted([y.split("/")[3].split(".")[0] for y in glob.glob("./../Dataset/*.pdb")])
print(len(filesList), filesList)

# all combinations
def find_specific_TSR_details(TSR): #correct()
    TSRname ='_'
    for t in TSR:
        TSRname = TSRname+t+'_'
    workbookname = subdirName+TSRname+'all_combinations_29_35.xlsx'
    writer = pd.ExcelWriter(workbookname, engine='xlsxwriter')
    
    avgfile = open(subdirName+TSRname+'all_combinations_average.txt', 'w')
    avgfile.write('prot\tavgTheta\tavgMaxDist\tstdTheta\tstdMaxDist\n')
    
    absent = []
    for i in filesList:    
        dfDetails= pd.DataFrame(columns=['key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'])
        c=0
        with open(subdirName+i+'.triplets_29_35', 'r') as file: 
            for r in file: # row is list
                if not r:
                    continue
                if r.split()[1]==TSR[0] and r.split()[3]==TSR[1] and r.split()[5]==TSR[2]: 
                     dfDetails.loc[c]= r.split()  # adding a row
                     c=c+1
                if r.split()[1]==TSR[0] and r.split()[3]==TSR[2] and r.split()[5]==TSR[1]: 
                     dfDetails.loc[c]= r.split()  # adding a row
                     c=c+1
                if r.split()[1]==TSR[1] and r.split()[3]==TSR[0] and r.split()[5]==TSR[2]: 
                     dfDetails.loc[c]= r.split()  # adding a row
                     c=c+1                    
                if r.split()[1]==TSR[1] and r.split()[3]==TSR[2] and r.split()[5]==TSR[0]: 
                     dfDetails.loc[c]= r.split()  # adding a row
                     c=c+1
                if r.split()[1]==TSR[2] and r.split()[3]==TSR[0] and r.split()[5]==TSR[1]: 
                     dfDetails.loc[c]= r.split()  # adding a row
                     c=c+1                    
                if r.split()[1]==TSR[2] and r.split()[3]==TSR[1] and r.split()[5]==TSR[0]: 
                     dfDetails.loc[c]= r.split()  # adding a row
                     c=c+1   
        print (i, c)
        if( c == 0):
            absent.append(i)
        dfDetails.theta = dfDetails.theta.astype(float)
        dfDetails.maxDist = dfDetails.maxDist.astype(float)
        c = c+1 # just to give one more line gap
        
        avgtheta = round(dfDetails['theta'].mean(),4)
        avgmaxdist = round(dfDetails['maxDist'].mean(),4)
        stdtheta = round(dfDetails['theta'].std(),4)
        stdmaxdist = round(dfDetails['maxDist'].std(),4)
        
        avgfile.write(i+'\t'+str(avgtheta)+'\t'+str(avgmaxdist)+'\t'+str(stdtheta)+'\t'+str(stdmaxdist)+'\n')
        
        dfDetails.loc[c,'theta']='avgTheta='+str(avgtheta)
        dfDetails.loc[c,'maxDist']='avgmaxDist='+str(avgmaxdist)
        c=c+2
        dfDetails.loc[c,'theta']='stdTheta='+str(stdtheta)
        dfDetails.loc[c,'maxDist']='stdmaxDist='+str(stdmaxdist)
        c=c+1
        dfDetails.to_excel(writer, sheet_name=i,index=False, startrow=0, startcol=0)
        avgfile.write(i+'\t'+str(avgtheta)+'\t'+str(avgmaxdist)+'\t'+str(stdtheta)+'\t'+str(stdmaxdist)+'\n')

    writer.save()
    avgfile.close()
    print("Proteins do not have this TSR: ",len(absent),absent)

def find_key_for_a_TSR(prot, TSR):
    TSRname ='_'
    for t in TSR:
        TSRname = TSRname+t+'_'
    f_out = open(subdirName+prot+'_'+TSRname+'_details.txt','w')
    file = subdirName+prot+'.triplets_29_35'
    with open(file,'r') as f:
        for lines in f:
            if not lines:
                continue
            line = lines.split()
            if line[2]==TSR[1] and line[4]==TSR[3] and line[6]==TSR[5]:
                print(line)
                f_out.write(lines)
            if line[2]==TSR[1] and line[4]==TSR[5] and line[6]==TSR[3]:
                print(line)
                f_out.write(lines)
            if line[2]==TSR[3] and line[4]==TSR[5] and line[6]==TSR[1]:
                print(line)
                f_out.write(lines)
            if line[2]==TSR[3] and line[4]==TSR[1] and line[6]==TSR[5]:
                print(line)
                f_out.write(lines)
            if line[2]==TSR[5] and line[4]==TSR[1] and line[6]==TSR[3]:
                print(line)
                f_out.write(lines)
            if line[2]==TSR[5] and line[4]==TSR[3] and line[6]==TSR[1]:
                print(line)
                f_out.write(lines)
    f_out.close()


def find_specific_key_details(protList, key): #correct()
    df_all = pd.DataFrame(columns=['key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2','prot'])
    outfile = subdirName+'_'+str(key)+'_details_29_35.txt'
    for p in protList:
        print(p)
        df = pd.read_csv(subdirName+p+'.triplets_29_35', sep='\t', names = ['key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'])
        df = df.loc[df['key'] == key]
        df['prot'] = p
        df_all= pd.concat([df_all,df])
        #print(df)
    print(df_all)
    df_all.to_csv(outfile, index=None, sep='\t')
    
def get_keyList():
    df = pd.read_csv('104Carboxypeptidases__HIS_GLU_HIS__details_specific_conditioning_29_35_maxdist8_new.txt', sep='\t', names = ['prot','key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'])
    keylist = (list(set(df['key'].tolist())))
    return keylist

def find_a_key_freq_in_each_protein(key, protList):
    print(key, len(protList))
    f_out = open(subdirName+str(key)+"_SER_sample_1_each_protein_freq.txt", "w")
    for p in sorted(protList):
        freq = 0
        with open(subdirName+p+".triplets_29_35", "r") as f:
              for lines in f:
                  if not lines.strip():
                      continue
                  line = lines.strip().split()
                  if len(line) != 20:
                      continue

                  if int(line[0].strip()) == key:
                      print(lines)
                      freq += 1
              f_out.writelines([str(key)+"\t"+p+"\t"+str(freq)+"\n"])
    f_out.close()
        
# func calls
#find_a_key_freq_in_each_protein(6362130, filesList)

find_specific_TSR_details(['ASP','HIS','SER']) #enter as list the TSR you want to find
find_specific_TSR_details(['ASP','HIS','THR']) #enter as list the TSR you want to find
find_specific_TSR_details(['GLU','HIS','SER']) #enter as list the TSR you want to find
find_specific_TSR_details(['GLU','HIS','THR']) #enter as list the TSR you want to find
find_specific_TSR_details(['ALA','VAL','LEU']) #enter as list the TSR you want to find
find_specific_TSR_details(['TRP','TYR','PHE']) #enter as list the TSR you want to find

#find_key_for_a_TSR('1ACB', ['HIS','57','SER','195','ASP','102'])

#find_specific_key_details(filesList, 6362130)
#print(get_keyList()) # get list of keys which has HIS_GLU_HIS with specific condition    

print ("code end")
  
