#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:29:08 2022

@author: schaferjw

Generate figures that show the change in the distribution of signal and noise as more subfamily alignments are incorporated
into the superimposition and more predictions are passed to the density based filter.

"""

import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from modules.DB_SORT import DB_SCAN
from modules.PLOT import PLOT
from modules.SUPER import Data
class Decomp():
    def __init__(self,files,pred_type):
        #load predictions from all subfamily alignements 
        self.cwd =  os.getcwd()
        result = [[0,j] for j in range(0,len(files)+1) if files[0:j]]
        if pred_type == 'FULL': 
            result = []
            for j in range(0,len(files)):
                if files[0:j+1]:
                    if 'gmn' in files[j] and 'msatr' in files[j+1]:
                        result.append([0,j+1])
                    if 'gmn' not in files[j-1] and 'msatr' in files[j]:
                        result.append([0,j])
        p_iter = {}
        for idx in result:
            p = Data(files,select=[pred_type,idx[0],idx[1]])
            p_iter[str(idx[1])] = p.Partial()
        self.p_iter = p_iter

    def Pred(self,msa,xcontact,msas,df,df_dist,pdb1,pdb2,Prediction_Offset,alignment,name,query_length):
        
        df_tot = df.append(df_dist)
        x_first =  df_tot[df_tot['Fold'] == pdb1]
        x_first =  x_first[['j','i']].to_numpy()
        x_second = df_tot[df_tot['Fold'] == pdb2]
        x_second = x_second[['i','j']].to_numpy()
        x_common = df_tot[df_tot['Fold'] == 'both']
        x_common = x_common[['j','i']].to_numpy()       
        
        gmn_count = open(f"{name}_super.txt","w")
        gmn_count.write("slice,percent,common,{:},{:},noise \n".format(pdb1,pdb2))
        gmn_count.close()
        
        m = sum(1 for line in open('{:}/{:}'.format(msa[:-4],msas[0])) if '>' not in line)
        output = [sum(1 for line in open('{:}/{:}'.format(msa[:-4],file)) if '>' not in line) for file in msas]
        output = ((7.5-1.5)/(np.max(output)-np.min(output)))*(output-np.min(output)) + 1.5
        keys = list(self.p_iter.keys())
        for key,msa_idx,N in zip(keys,msas,list(reversed(output))):
            n = sum(1 for line in open('{:}/{:}'.format(msa[:-4],msa_idx)) if '>' not in line)
            df_super = self.Most_Probable(self.p_iter[key],N_r=N)
            if Prediction_Offset != 'n':
                df_super['i'] = df_super['i'] + int(Prediction_Offset)
                df_super['j'] = df_super['j'] + int(Prediction_Offset)
            
            db = DB_SCAN(df_super,xcontact,alignment)
            alignment = db.Return_Alignment()
            df_sorted = db.Run(x_first,x_second,x_common,msa,query_length,'y')
            df_sorted.loc[df_sorted.sort == 'pdb_1', 'sort'] = pdb1
            df_sorted.loc[df_sorted.sort == 'pdb_2', 'sort'] = pdb2
            counts = df_sorted['sort'].value_counts()
            if pdb1  not in counts: counts[pdb1]  = 0
            if pdb2  not in counts: counts[pdb2]  = 0
            if 'common' not in counts: counts['common'] = 0
            if 'noise'  not in counts: counts['noise']  = 0
            gmn_count = open(f"{name}_super.txt", "a")  # append mode
            gmn_count.write("{:},{:},{:},{:},{:},{:} \n".format(key,float(1-(n/m)),counts['common'],counts[pdb1],counts[pdb2],counts['noise']))
            gmn_count.close()

        df_full  = pd.read_csv(f'{name}_super.txt')
        df_full.plot(x='percent',y=[pdb1,pdb2,'common','noise '] ,style='.-')
        plt.savefig(f'{name}_super.png')
           
    def Most_Probable(self,mtx,N_r=2):
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
