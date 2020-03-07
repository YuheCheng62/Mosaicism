import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.linalg import null_space
data=pd.read_csv('20191028_Table_for_analysis.csv',index_col=0)
Nanmap=data.isna()#change NaN to 0
for row in data.index:
    for column in data.columns:
        if Nanmap.loc[row,column]==True:
            data.loc[row,column]=0
AF=np.array(data.values)[:,0:11]#Taking the first 11 snp site
x=null_space(AF)#Find null space of Allele frequency matrix as the normal vector of hyperplane
snp_list=[]
print(np.matmul(AF,x))#for each point belongs to {0,1}^n, calculate the dot product to check if it is on the hyperplane 
for ii in range(2047):
    bin_str='{0:b}'.format(ii)
    point=np.zeros(11)
    for jj,char in enumerate(bin_str):
        if char=='1':
            point[11-len(bin_str)+jj]=1
    dot_pro=np.dot(np.array(point),x)
    if dot_pro<0.002 and dot_pro>-0.002:
        print(dot_pro)
        print(point)
        snp_list.append(point)
# Combine all binary points that are on the hyperplane
snp_matrix=np.array(snp_list)


