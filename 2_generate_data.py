"""
Created on Mon Sep 30 21:08:26 2019

@author: Titli

Decide central, lower, higher point and add d1, d2 in files | folder /d1d2 | file - '_distance_added.csv'
"""
import glob, os, time
import pandas as pd 
import argparse
import multiprocessing
from joblib import Parallel, delayed
import numpy as np 

start_time =time.time()
#python code.py input_dir inputfile_extension    
parser = argparse.ArgumentParser()
parser.add_argument("data_dir", default='./../Dataset/lexiographic/', help="enter path to the input file", type=str) # ./../Dataset/Intra_data_negetive_key/
parser.add_argument("extension", default='.triplets_29_35', help="enter file extension") 
args = parser.parse_args()

data_dir = args.data_dir # input files location
subdir = data_dir+'d1d2/'
if not os.path.exists(subdir):
    os.makedirs(subdir)

data_location = len(data_dir.split("/"))-1
print(data_dir, "\n", data_dir.split("/"))
print(data_location)

files = [x.split("/")[data_location].split(".")[0] for x in glob.glob(data_dir+"*"+args.extension)]
print(len(files), files)

# helper function
def calculate_distance(pos0, pos1, pos2):
    # determine central, lower and higher points
    positions = sorted([int(pos0), int(pos1), int(pos2)])
    central_point = positions[1] #find the central point from positions
    lower_point = positions[0]
    higher_point = positions[2]
    #print(positions, central_point, lower_point, higher_point)
    
    # calculate d1, d2 from cental point       
    d1 = abs(central_point-lower_point) # distance between central point and lower position point
    d2 = abs(central_point-higher_point) # distance between central point and lower position point
    return d1,d2

def prepare_distance_data(file):
    print(file)
    start = time.time() 
    # get the freq of each TSR and only consider the TSRs having freq >= 2
    #df = pd.read_csv(data_dir+file+args.extension, sep="\t",usecols=[0,1,2,3,4,5,6], names=['key','amino0','pos0','amino1','pos1','amino2','pos2'])
    df = pd.read_csv(data_dir+file+args.extension, sep="\t",header=0, usecols=[0,1,2,3,4,5,6])

    #counts = df['key'].value_counts()
    #df1 = df[df['key'].isin(counts[counts >= 2].index)] # get all keys(including all coourences of each key) with freq>=2 
    #print(df1)
    # for all key occurence, calculate (d1, d2)
    d1_list = []
    d2_list = []
    for index,row in df.iterrows(): #iterate over df rows with iterrows()
        #if (row['pos0']=='pos0' or row['pos1']=='pos1' or row['pos2']=='pos2'):
            #continue
        #print(row['pos0'], type(row['pos0']))
        #print(row['pos1'], type(row['pos1']))
        #print(row['pos2'], type(row['pos1']))

        d1,d2 = calculate_distance(row['pos0'], row['pos1'], row['pos2'])
        d1_list.append(d1)
        d2_list.append(d2)

    #print(len(df1), len(d1_list), len(d2_list))
    df['d1'] = d1_list
    df['d2'] = d2_list
    df['protein'] = file
    df.to_csv(subdir+file+'_distance_added.csv',index=None, sep='\t', mode='w')
    print(file," Distance Data creation time(sec): ", round(time.time()-start, 4))

num_cores = multiprocessing.cpu_count()
print("#of cores = ", num_cores)
Parallel(n_jobs=num_cores, verbose=50)(delayed(prepare_distance_data)(fileName)for fileName in files)

print("CODE END\t TIME=(min)", np.round((time.time()-start_time)/60,2))