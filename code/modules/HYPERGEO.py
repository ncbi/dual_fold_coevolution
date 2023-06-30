#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 10:03:58 2022

@author: schaferjw

run hypergeometric test to calculate a P-value corresponding to the statistical significance of the predictions
compared to the original coevolutionary calculations (GREMLIN/MSATransformer)
"""

from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

class HYPERGEO():
    def __init__(self,pdb1,pdb2,msa,df_gmn_msatr):
	#load information for calculating p-value
        df_xcontact = pd.read_csv('{:}/df_xcontact.csv'.format(msa[:-4]),index_col=0)
        df_xcontact_dist = pd.read_csv('{:}/df_dist_xcontact.csv'.format(msa[:-4]),index_col=0)
        xcontact = pd.concat([df_xcontact,df_xcontact_dist])
        pairs = xcontact[['i','j']].to_numpy()
        points = np.array([(0,1),(1,0),(0,-1),(-1,0),(1,1),(-1,-1),(-1,1),(1,-1)] )
        
        for pair in pairs:
            off = np.concatenate([pair+points for pair in pairs],axis=0)
        new_array = [tuple(row) for row in off]
        uniques = np.unique(new_array,axis=0)
        self.K = len(uniques)     
        
        self.df_gmn_msatr = df_gmn_msatr
        self.N = len(df_gmn_msatr.index)

    def P_value(self, df_super, df_original, df_orig_msatr, N, x_first, x_second, x_common):
        # calculat p-value using hypergeometeric test 
        k=len(df_super[(df_super['group'] != -1) & (df_super['sort'] != 'noise')]) - len(self.df_gmn_msatr[(self.df_gmn_msatr['group'] != -1) & (self.df_gmn_msatr['sort'] != 'noise')])
        n=len(df_super[df_super['group'] != -1].index) - len(self.df_gmn_msatr[self.df_gmn_msatr['group'] != -1].index)

        N=N**2 - self.N
        K=self.K
        
        rv = hypergeom(N, K, n)
        x = np.arange(0, n+1)
        pmf_pred = rv.pmf(x)
        
        print('P-value: {:}'.format(sum(pmf_pred[k:])))