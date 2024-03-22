# -*- coding: utf-8 -*-

"""
Created on 20/02/2023

@author: terence
"""

import numpy as np
import pandas as pd
import numpy.matlib
import os
### arguments
#finname=sys.argv[1]
#sN=sys.argv[2]

### constants


dir = "C:\\Users\\legrandt\\Nextcloud2\\PostDoc\\Postdoc_CCMAR\\Projet\\Polynesia_oyster\\DATA\\"
site_dir = "Genetics\\"
mat_dir = "Connectivity_matrices_processed\\filtered\\"
mat_save = "Connectivity_matrices_processed\\temp_mean\\"

# Load site data
site_data = pd.read_csv(dir + site_dir + "df_islands.csv")

# Get all files names
list_mat = os.listdir(dir + mat_dir)

N=len(site_data)

for x in range(0,len(list_mat)):
    print("--> "+str(round(x/len(list_mat),3))+"%")
                                                                                                                                            
    matrix=np.zeros((N,N))                                                                                                                                                                                                  
### reading matrix  

    p=open(dir + mat_dir + "\\" + list_mat[x])  
    #p=open('test.matrix', "r") #test                                                                                                                                           
# =============================================================================
    
    for line in p:                                                                                                                       
        matrix[int(line.split()[0]),int(line.split()[1])] = float(line.split()[2])                     

### instrength calculation
        
    instr=np.sum(matrix, axis=0)

### istrength matrix size
    
    instr=np.matlib.repmat(instr,N,1)
    
### normalizing matrix
       
    matrix=np.divide(matrix,instr)
    
    ## Replace NaN with zero

    matrix=np.nan_to_num(matrix,nan=0.0)
    
#    for i in range (0,N):
#        for j in range (0,N):
#            if instr[i] > 0. :
#                matrix[i,j] = matrix[j,i]/instr[i]
#    
    p.close()

#### writing matrix
    d=open(dir + mat_save + list_mat[x],'w')
    for i in range (0,N):
            for j in range (0,N):
                if matrix[i,j] > 0. :
                    d.write(str(i) + ' ' + str(j) + ' ' + str(matrix[i,j]) + '\n')                                                                                                                                                       
    d.close()
    
# =============================================================================