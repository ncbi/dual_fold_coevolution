#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:29:08 2022

@author: schaferjw

Functions needed to create the superimposition of coevolutionary predictions

"""

import numpy as np
import pandas as pd
import os
from os import listdir
from os.path import isfile, join

class Data():
    def __init__(self,files,select=['NA']):
        #load predictions from all subfamily alignements
        self.cwd =  os.getcwd()

        if select[0] != 'NA':files = files[select[1]:select[2]]
        
        coevolution = [self.Load_Numpy(file) for file in files]
        coevolution = [self.Normalize(mtx) for mtx in coevolution]
        coev_super = self.Superposition(coevolution)
        
        #self.mtx predictions will be sent through the rest of the pipeline
        self.mtx = coev_super
        
        if select[0] != 'NA':self.partial = coev_super
        
    def Load_Numpy(self,file):
        zero = np.loadtxt(file,delimiter=",", dtype=float)
        return zero
    
    def Superposition(self,mtxs):
        mtxs = np.array(mtxs)
        mtx_sup = np.zeros((mtxs.shape[1],mtxs.shape[2]))
        for i in range(mtxs.shape[1]):
            for j in range(mtxs.shape[2]):
                mtx_sup[i,j] = np.average(mtxs[:,i,j])
        return mtx_sup
                
    def Normalize(self,mtx):
        i,j,z = [],[],[]
        for a in range(mtx.shape[0]):
            for b in range(mtx.shape[1]):
                i.append(a)
                j.append(b)
                z.append(mtx[a,b])
        z = np.array(z)

        output = ((1.0-0.0)/(np.max(z)-np.min(z)))*(z-np.min(z)) + 0.0

        zero = np.zeros(mtx.shape)
        for idx in range(len(output)):
            zero[i[idx],j[idx]] = output[idx] 
        return zero
    
    def Most_Probable(self,N_r=2):
        mtx=self.mtx
        #N_r dictates the number of high confidence points that will be returned
        N_r = int(N_r*mtx.shape[0])

        N = int(mtx.shape[0])
        i,j,z = [],[],[]
        for a in range(mtx.shape[0]):
            for b in range(mtx.shape[1]):
                i.append(a)
                j.append(b)
                z.append(mtx[a,b])
        df = pd.DataFrame({'i': i, 'j': j,'zscore': z})
        top = df.loc[df['j'] - df['i'] >= 3].sort_values("zscore",ascending=False)
        temp = top.head(N_r)

        return temp
    
    def Partial(self):return self.partial