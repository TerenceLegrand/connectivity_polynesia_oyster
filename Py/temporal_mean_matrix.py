"""
Created on 20/02/2023

@author: terence
"""

#%reset

import numpy as np
import pandas as pd
import numpy.matlib
import os

dir = "C:\\Users\\legrandt\\Nextcloud2\\PostDoc\\Postdoc_CCMAR\\Projet\\Polynesia_oyster\\DATA\\"
mat_dir = "Connectivity_matrices_processed\\back_norm\\"
mat_save = "Connectivity_matrices_processed\\temp_mean\\"

# Load site data
site_data = pd.read_csv(dir + "df_islands.csv")

# Get all files names
list_mat = os.listdir(dir + mat_dir)

N=len(site_data)

matrix=np.zeros((N,N)) 

# Add all matrices
for x in range(1,len(list_mat)): # from 1 because list_mat[1] is a folder
    print("--> "+str(round(x/len(list_mat),3))+"%")
                                                                                                                                            
                                                                                                                                                                                                     

    p=open(dir + mat_dir + "\\" + list_mat[x])  

    for line in p:
        matrix[int(line.split()[0]),int(line.split()[1])] = matrix[int(line.split()[0]),int(line.split()[1])] + float(line.split()[2])

    p.close()

# Temporal normalisation
matrix=np.true_divide(matrix,len(list_mat))

# Backward normalisation 
instr=np.sum(matrix, axis=0)
instr=np.matlib.repmat(instr,N,1)
matrix=np.divide(matrix,instr)
matrix=np.nan_to_num(matrix,nan=0.0)

# Transpose matrix for multi-generation connectivity
matrix=matrix.transpose()

#### writing matrix
d=open(dir + mat_save + "back_temporal_mean_2010_2019" + ".matrix",'w')
for i in range (0,N):
        for j in range (0,N):
            if matrix[i,j] > 0. :
                d.write(str(i) + ' ' + str(j) + ' ' + str(matrix[i,j]) + '\n')                                                                                                                                                       
d.close()

