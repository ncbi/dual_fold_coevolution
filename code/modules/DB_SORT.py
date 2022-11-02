#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:32:20 2022

@author: schaferjw
"""

import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
from kneed import KneeLocator
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from os import listdir
from os.path import isfile, join

class DB_SCAN():
    def __init__(self,df,xcontact,alignment='n'):
        #search for optimal alignment to dual fold xcrytal contacts
        if alignment == 'n':
            result = []
            for o in range(-75,76):
                temp = df.copy()
                temp['i'] = temp['i'] + o
                temp['j'] = temp['j'] + o
                
                count = 0
                for index,row in temp.iterrows():
                    if (int(row['i']),int(row['j'])) in xcontact or (int(row['j']),int(row['i'])) in xcontact:
                        count+=1
                    else:flag=False
                if flag == False:result.append((o,count))
            best = result[0]
            temp = df.copy()
            for i in result:
                if i[1] > best[1]:
                    best = i
            temp['i'] = temp['i'] + best[0]
            temp['j'] = temp['j'] + best[0]

            self.align = temp
            self.alignment = best[0]
        else:
            self.alignment = int(alignment)
            df['i'] = df['i'] + int(alignment)
            df['j'] = df['j'] + int(alignment)
            self.align = df

    def Return_Alignment(self): return self.alignment
    
    def DB(self,contacts,knee,min_samples):
        #use DBscan algorithm on superposition to eliminate noise
        db = DBSCAN(eps=knee, min_samples=min_samples).fit(contacts)
        core_mask = np.zeros_like(db.labels_,dtype=bool)
        core_mask[db.core_sample_indices_]=True
        
        points = {}
        for label in set(db.labels_):
            class_mask = db.labels_ == label
            clusters = contacts[class_mask & core_mask]
            clusters2 = contacts[class_mask & ~core_mask]
            temp = [[contact[0],contact[1]] for contact in clusters]
            if np.any(temp) == True: points[label] = temp
            temp = [[contact[0],contact[1]] for contact in clusters2]
            if np.any(temp) == True: points[label] = temp
        
        return points
    
    def Run(self,x_first,x_second,x_common,msa, no_filter='n'):
        
        #define input to find optimal parameters 
        contacts = self.align[['i','j']].to_numpy()
        knee = self.DB_opt(contacts)
        num = 50
        param,best = [],[]
	#search around the "knee" of the graph to construct an ROC curve
        for a in np.linspace(knee*0.25,knee*1.1,num=num):
            points_sort = {}
            min_samples = 2 
            points = self.DB(contacts,a,min_samples)
            for key in points:
                if len(points[key]) > 2 and no_filter=='n': #remove sparse clusters of 1 or 2 points 
                    points_sort[key] = self.GMN_decomp(np.array(points[key]), x_first, x_second, x_common)
                elif no_filter == 'y':
                    points_sort[key] = self.GMN_decomp(np.array(points[key]), x_first, x_second, x_common)
            df_sorted = pd.DataFrame(columns=["i","j","apc","zscore","gmn_x_freq","group","sort","r"])
            for key in points_sort:
                for contact in points_sort[key]:
                    temp = self.align.loc[(self.align['i'] == contact[0]) & (self.align['j'] == contact[1])]
                    temp['gmn_x_freq'] = temp['zscore'] * temp['freq']
                    temp['group'] = key
                    temp['sort']  = contact[2]
                    temp['r']     = contact[3]
                    df_sorted = df_sorted.append(temp)
           
	    #find tpr and fpr for this iteration 
            if no_filter == 'n':
                p = df_sorted.drop(df_sorted[df_sorted.group == -1].index).copy()
                n = df_sorted.drop(df_sorted[df_sorted.group != -1].index).copy()
                if p.empty != True:
                    if 'noise' in p['sort'].unique():tp = len(p.index)-p['sort'].value_counts().loc[['noise']].item() #total-noise
                    else:tp=len(p.index)
                    if 'noise' in n['sort'].unique():fn = len(n.index)-n['sort'].value_counts().loc[['noise']].item() #total-noise
                    else:fn=len(n.index)
                    tpr = tp/(tp+fn)
                    if 'noise' in p['sort'].unique():fp = p['sort'].value_counts().loc[['noise']].item() #noise
                    else:fp=0
                    if 'noise' in n['sort'].unique():tn = n['sort'].value_counts().loc[['noise']].item() #noise
                    else:tn=0
                    fpr = fp/(fp+tn)
                
                    if not param:param = [[fpr,tpr,a,min_samples,len(df_sorted[df_sorted['group'] != -1].index)]]
                    if param[-1][0] != fpr and param[-1][1] != tpr:
                        param.append([fpr,tpr,a,min_samples,len(df_sorted[df_sorted['group'] != -1].index)])
            else:param = [[1,1,a,min_samples,len(df_sorted[df_sorted['group'] != -1].index)]]

        files = [i for i in listdir('{:}'.format(msa[:-4])) if isfile(join('{:}'.format(msa[:-4]), i))]
        files = [file for file in files if file[0:2] == 'df']
        files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])
        df_orig = pd.read_csv('{:}/{:}'.format(msa[:-4],files[0]),index_col=0)
        min_pred = len(df_orig.index)
        
        for i in range(len(param)-1): #define best point on ROC curve
            y = param[i+1][1] - param[i][1]
            x = param[i+1][0] - param[i][0]
            if y/x < 1 and param[i][4] > min_pred:
                best = param[i]
                break
        
        if not best: # if there is no improvement return the last point on the ROC curve (retain the most information)
            best = param[-1]
        
        if no_filter == 'n': #create png of ROC curve
            fig, axs = plt.subplots(figsize=(6, 6))
            df_param = pd.DataFrame(param,columns=['fpr','tpr','eps','min_samples','min_pred'])
            df_param.plot.scatter(x='fpr',y='tpr', ax=axs)
            df_best = pd.DataFrame([best],columns=['fpr','tpr','eps','min_samples','min_pred'])
            df_best.plot.scatter( x='fpr',y='tpr', ax=axs, color='r')
            axs.set_ylim(0,1)
            axs.set_xlim(0,1)
            fig.savefig("roc.png")
            
        self.eps = best[2:]
        points_sort = {}
        points = self.DB(contacts,best[2],best[3]) #use optimal parameters to get final set of filtered predictions
        for key in points:
            if len(points[key]) > 2 and no_filter=='n': 
                points_sort[key] = self.GMN_decomp(np.array(points[key]), x_first, x_second, x_common)
            elif no_filter == 'y':
                points_sort[key] = self.GMN_decomp(np.array(points[key]), x_first, x_second, x_common)
                
        df_sorted = pd.DataFrame(columns=["i","j","apc","zscore","gmn_x_freq","group","sort","r"])
        for key in points_sort:
            for contact in points_sort[key]:
                temp = self.align.loc[(self.align['i'] == contact[0]) & (self.align['j'] == contact[1])]
                temp['gmn_x_freq'] = temp['zscore'] * temp['freq']
                temp['group'] = key
                temp['sort']  = contact[2]
                temp['r']     = contact[3]
                df_sorted = df_sorted.append(temp)  
       
	#remove sequence predicitons that are outside of the xcrystal region 
        df_sorted = self.Trim_Lower(df_sorted)
        df_sorted = self.Trim_Upper(df_sorted)
        return df_sorted
    
    def Trim_Lower(self,df):
        for i in range(df['j'].astype('int32').min(),df['j'].astype('int32').max()):
            b=2*i
            df['temp'] = df['i'] + df['j'] - b
            check = df.loc[df['temp'] < 0]['sort'].value_counts().to_dict()
            if 'noise' not in check.keys() and bool(check) == True:break
            if 'noise' in check.keys() and len(check) > 1:break
        b=b-2
        df['temp'] = df['i'] + df['j'] - b
        df = df.drop(df[df.temp < 0].index)
        
        return df.drop(['temp'], axis=1)
    def Trim_Upper(self,df):
        for i in reversed(range(df['j'].astype('int32').min(),df['j'].astype('int32').max())):
            b=2*i
            df['temp'] = df['i'] + df['j'] - b
            check = df.loc[df['temp'] > 0]['sort'].value_counts().to_dict()
            if 'noise' not in check.keys() and bool(check) == True:break
            if 'noise' in check.keys() and len(check) > 1:break
        b=b+2
        df['temp'] = df['i'] + df['j'] - b
        df = df.drop(df[df.temp > 0].index)
        return df.drop(['temp'], axis=1)

    def EPS(self):return self.eps
    def DB_opt(self,predictions):
        #Algorithm for finding the "elbow" of a graph, this gives a reasonable starting point
	#original work: https://ieeexplore.ieee.org/document/5961514
        nbrs = NearestNeighbors(n_neighbors=4)
        nbrs_fit = nbrs.fit(predictions)
        distances, indices = nbrs_fit.kneighbors(predictions)
        distances = np.sort(distances, axis=0)
        distances = distances[:,1]
        
        x = [i for i in range(len(distances))]
        y = distances
        kneedle = KneeLocator(x, y, S=8, curve="convex", direction="increasing")
        
        return kneedle.knee_y

    def Sort_fold(self,pdb_1_r,pdb_2_r,common_r):
        #define which xcrystal structure contact a predicion matches (note: predictions are allowed to be off by +-1)
        if pdb_1_r[0] > 1.5 and pdb_2_r[0] > 1.5 and common_r[0] > 1.5:return ['noise', min(pdb_1_r[0],pdb_2_r[0],common_r[0])]
        else:
            if common_r[0] == 0:                                      return ['common', min(pdb_1_r[0],pdb_2_r[0],common_r[0])]
            if common_r[0]< pdb_1_r[0] and common_r[0] < pdb_2_r[0]:  return ['common', min(pdb_1_r[0],pdb_2_r[0],common_r[0])]
            if pdb_1_r[0]== pdb_2_r[0] and pdb_1_r[0]== common_r[0]:  return ['common', min(pdb_1_r[0],pdb_2_r[0],common_r[0])]
            if pdb_1_r[0]== pdb_2_r[0]:                               return ['common', min(pdb_1_r[0],pdb_2_r[0],common_r[0])]
            if pdb_1_r[0] < pdb_2_r[0] and pdb_1_r[0] <= common_r[0]: return ['pdb_1',  min(pdb_1_r[0],pdb_2_r[0],common_r[0])]
            if pdb_2_r[0] < pdb_1_r[0] and pdb_2_r[0] <= common_r[0]: return ['pdb_2',  min(pdb_1_r[0],pdb_2_r[0],common_r[0])]
    
    def GMN_decomp(self,group,x_first,x_second,x_common):
	#define the distance between predicted contacts and closest xcrystal contact
        Tree_first  = scipy.spatial.cKDTree(x_first)
        Tree_second = scipy.spatial.cKDTree(x_second)
        Tree_common = scipy.spatial.cKDTree(x_common)
        
        first_iter  = [Tree_first.query(group[idx].reshape((1,2)), k=1,distance_upper_bound=200) for idx in range(len(group))]
        second_iter = [Tree_second.query(group[idx].reshape((1,2)),k=1,distance_upper_bound=200) for idx in range(len(group))]
        common_iter = [Tree_common.query(group[idx].reshape((1,2)),k=1,distance_upper_bound=200) for idx in range(len(group))]
        sort_iter   = [self.Sort_fold(idx[0],idx[1],idx[2]) for idx in zip(first_iter,second_iter,common_iter)]        
        
        return [[idx[0][0],idx[0][1],idx[1][0],idx[1][1]] for idx in zip(group,sort_iter)]
