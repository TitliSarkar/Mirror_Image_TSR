# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:48:20 2021

@author: titli
"""
import os, glob, time
import pandas as pd
import numpy as np 

pd.set_option('display.max_columns', None)

start = time.time()

file = 'C:\\D\\CODES\\Python_Code\\SET1\\Sequence_to_TSR\\result_key_percentage_lex_mirror_maxdist15_20_0.csv'
 
print(file)
df = pd.read_csv(file, header=0, sep="\t")
print(df.columns)
#df = df.sort_values(by='percentage(total_unq/total_dataset)', ascending=False)

maxrow = df.loc[[df['percentage(total_unq/total_dataset)'].argmax()]]
print(maxrow)

print("Time (min) =", np.round((time.time()-start)/60, 2))