import numpy as np
import pandas as pd
import os.path

dataset = os.getcwd().split("/")[-2]
print(dataset)

dir_ = "./../Dataset/lexiographic/"

def calculate_ratio(infile, outfile):
    df= pd.read_csv(infile,header=0, sep='\t')
    #print(df)
    f_out = open(outfile, 'w')

    total = np.sum(df['total'].tolist())
    total_u_u_u = np.sum(df['total u_u_u'].tolist())
    total_u_u_c = np.sum(df['total u_u_c'].tolist())
    total_u_c_c = np.sum(df['total u_c_c'].tolist())
    total_c_c_c = np.sum(df['total c_c_c'].tolist())

    N = total
    N1 = total_u_u_u
    N2 = total_u_u_c
    N3 = total_u_c_c
    N4 = total_c_c_c

    total_singleton_aa_percentage = np.round((N1/N)*100, 2)
    total_grouped_aa_percentage = np.round(((N2+N3+N4)/N)*100, 2)

    ratio1 = np.round(((N2 + N3 + N4)/ N1),2)
    ratio2 = np.round(((N4+ N3)/N1),2)
    ratio3 = np.round(((N4 + N3)/(N1+N2)),2)
    ratio4 = np.round((N4/N1),2)
    ratio5 = np.round((N4/ (N1 + N2 + N3)),2)

    tmp = [total,total_u_u_u,total_u_u_c,total_u_c_c,total_c_c_c, ratio1, ratio2, ratio3, ratio4, ratio5, total_singleton_aa_percentage, total_grouped_aa_percentage]
    line = "\t".join(str(x) for x in tmp)
    
    f_out.write('total'+'\t'+'total_u_u_u'+'\t'+'total_u_u_c'+'\t'+'total_u_c_c'+'\t'+'total_c_c_c'+'\t'+'ratio1'+'\t'+'ratio2'+'\t'+'ratio3'+'\t'+'ratio4'+'\t'+'ratio5'+'\t'+'total_singleton_aa_percentage'+'\t'+'total_grouped_aa_percentage'+'\n')
    
    f_out.write(line)
    f_out.close()

calculate_ratio(dir_+dataset+'_TSR_percentage.txt',dir_+dataset+'_ratios_over_collection.txt')
calculate_ratio(dir_+dataset+'_TSR_percentage_commonKeys.txt',dir_+dataset+'_ratios_over_collection_commonKeys.txt')


