# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 12:16:16 2019

@author: Titli
""" 
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

data_dir = '.\\..\\Dataset\\'

df = pd.read_csv(data_dir+'ssim_test_old.image.txt', header=0, sep="\t")
#df = df.drop('prot', axis=1)
df = df.iloc[:, :-1]
indices = df['prot'].tolist()
df = df.set_index('prot')
#df = df.drop('prot', axis=1)
print(df, df.shape)

# dendrogram
Z = linkage(df, 'average')
plt.title('dendrogram_ssim_test_old.image')
dendrogram(Z, labels=indices, orientation='left')
plt.savefig(data_dir+'dendrogram_ssim_test_old.image.jpg', bbox_inches='tight')
plt.clf()

# dendrogram with heatmap
sns.set(font_scale=0.8)
sns_plot = sns.clustermap(df, method="average",xticklabels=True, yticklabels=True)#, figsize=(50,50))
sns_plot.savefig(data_dir+'heatmap_ssim_test_old.image.jpg')


