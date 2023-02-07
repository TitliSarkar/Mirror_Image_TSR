# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 13:46:43 2018

@author: Titli
"""

#from scipy.cluster.hierarchy import linkage

# Libraries

import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
import numpy as np
import pandas as pd


dataset = os.getcwd().split("/")[-2]
print(dataset)

dir_ = "./../Dataset/lexiographic/"

X = np.loadtxt(dir_+'JaccardSimilarity_'+dataset+'_29_35.txt')
Z = hierarchy.linkage(X, 'ward')

names = ['5DWB_RE', '1SX8_RE', '1BX6_PKABC', '1MRV_PKABC', '3NX8_Kinase', '4UY9_Kinase', '1H1B_Protease', '1HXF_Protease', '1SHL_Caspase', '4EHA_Caspase']
print(len(names), names)

plt.title('44PKABC')
plt.xlabel('sample index')
plt.ylabel('distance (Ward)')
hierarchy.dendrogram(Z, orientation="left", leaf_rotation=0, leaf_font_size=8, labels=names)
plt.savefig('plt1.png', dpi=320, format='png', bbox_inches='tight')

#print(len(names))
#df= pd.read_excel('Protein ID for TSR Manuscript 8-24-2018.xlsx',sheet_name='PKABC')
#df['prot_group'] = df['PDB IDs'] +'_'+ df['Name']
#names = df['prot_group'].tolist()
#names = ['5DWB_RE', '1SX8_RE', '1BX6_PKABC', '1MRV_PKABC', '3NX8_Kinase', '4UY9_Kinase', '1H1B_Protease', '1HXF_Protease', '1SHL_Caspase', '4EHA_Caspase']
#print(len(names), names)

#figure = ff.create_dendrogram(X, labels=names,orientation='left', linkagefun=lambda x: linkage(X, 'ward', metric='euclidean'))
#figure['layout'].update({'width':800, 'height':800})

#py.plot(figure, filename='dendrogram_aa_mix#1_29_35_aa_grp_0.jpeg')
#plotly.offline.iplot(figure, filename='dendrogram_aa_mix#1_29_35_offline.jpeg')