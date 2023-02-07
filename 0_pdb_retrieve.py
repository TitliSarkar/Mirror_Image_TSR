# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 12:39:15 2018

@author: Titli
"""
'''# Dr. Xu: open anaconda-prompt, then type
     pip install urllib'''
 
import urllib.request      
import multiprocessing
from joblib import Parallel, delayed
import pandas as pd 
import os, argparse

parser = argparse.ArgumentParser(description='Downloading pdb files.')
parser.add_argument('sample_Details', default='sample_details.csv',type=str, help='Enter input file for protein list to be downloaded. File format should be .csv file with first column having list of proteins with header=protein')
args = parser.parse_args()

out_dir = "./../Dataset/" # output dir
os.makedirs(out_dir, exist_ok=True)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

df = pd.read_csv(args.sample_Details)
protList = df['protein']
print(len(protList), protList)


def parallelcode(fname):
    print (fname)
    urllib.request.urlretrieve('http://files.rcsb.org/download/'+fname+'.pdb', out_dir+fname+'.pdb')

def getInputs():
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores, verbose=50)(delayed(parallelcode)(fname)for fname in protList)

getInputs()

print ("code end.")
