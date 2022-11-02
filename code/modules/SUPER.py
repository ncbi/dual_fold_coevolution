#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:29:08 2022

@author: schaferjw
"""

import numpy as np
import pandas as pd
import os
from os import listdir
from os.path import isfile, join

class GMN_Data():
    def __init__(self,file_in):
        #load predictions from all subfamily alignements into a dataframe
        self.cwd =  os.getcwd()
        self.file_in = file_in
        
        files = [i for i in listdir('{:}/{:}'.format(self.cwd,self.file_in[:-4])) if isfile(join('{:}/{:}'.format(self.cwd,self.file_in[:-4]), i))]
        raw_predictions = pd.DataFrame(columns=["i","j","apc","zscore","df"])
        count = 0
        for file in files:
            if file[0:6] == 'df_gmn' or file[0:8] == 'df_msatr' and file != 'df_super.csv':
                count+=1
                if raw_predictions.empty == True: 
                    raw_predictions = pd.read_csv('{:}/{:}'.format(self.file_in[:-4],file),index_col=0)
                    raw_predictions['df'] = file
                else:
                    temp = pd.read_csv('{:}/{:}'.format(self.file_in[:-4],file),index_col=0)
                    temp['df'] = file
                    raw_predictions = pd.concat([raw_predictions,temp],ignore_index=True)
        
        self.raw_predictions = raw_predictions
        self.count = count
    
    def Super(self):
	#sort through all predictions
        temp = self.raw_predictions.copy()
        dict_av = {}
        for index,row in temp.iterrows():
            av = pd.DataFrame(columns=["i","j","apc","zscore"])
            comp = temp.loc[[index]]
            if '{:}-{:}'.format(comp['i'][index],comp['j'][index]) not in dict_av:
                for index_2,row_2 in temp.iterrows():
                    if row_2['i'] == comp['i'][index] and row_2['j'] == comp['j'][index]:
                        app = temp.loc[[index_2]]
                        av = pd.concat([av,app])    
                        dict_av['{:}-{:}'.format(comp['i'][index],comp['j'][index])] = av

        pd.options.mode.chained_assignment = None  # default='warn'                        
        final = pd.DataFrame(columns=["i","j","apc","zscore"])
        for key in dict_av:
            position = dict_av[key]
            qid_slice = position['df']
            qid_slice = [row['df'] for index,row in position.iterrows()]
            qid_slice = sorted(qid_slice, key = lambda x: x.rsplit('.', 1)[0])
            if 'df_super.csv' in qid_slice:
                qid_slice.remove('df_super.csv')
            qid_slice = '-'.join(qid_slice)
            temp = position.iloc[[0]]
            apc = temp['apc']
            zscore = temp['zscore']

            for index,row in position.iterrows():
                apc = (row['apc'] + apc) 
                zscore = (row['zscore'] + zscore)
            apc = apc / len(position.index)
            zscore = zscore / len(position.index)
            temp['apc'] = apc
            temp['zscore'] = zscore
            temp['freq'] = position.shape[0] / self.count
            temp['df'] = qid_slice
            final = pd.concat([final,temp])

        path = self.cwd+'/{:}/df_super.csv'.format(self.file_in[:-4])
        final.to_csv(r'{:}'.format(path))
