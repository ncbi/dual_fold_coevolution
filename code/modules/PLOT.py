#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 15:58:10 2022

@author: schaferjw
"""
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class PLOT():
    def FS_plot(self,pdb1,pdb2,df_dist,df,msa,df_sorted):
        colors = {'both':'#767676'}
        counts = df_sorted['sort'].value_counts()
        if len(counts) == 4 and counts['pdb_1'] > counts['pdb_2']:
            df_sorted = df_sorted.rename({'i':'j', 'j':'i'}, axis=1)
            df_dist = df_dist.rename({'i':'j', 'j':'i'}, axis=1)
            df = df.rename({'i':'j', 'j':'i'}, axis=1)
            #Make sure the both is always gray
            colors = {'both':'#767676', '{:}'.format(pdb1):'#d8d8d8', '{:}'.format(pdb2):'#000000'}
        if len(counts) == 4 and counts['pdb_1'] < counts['pdb_2']:
            colors = {'both':'#767676', '{:}'.format(pdb2):'#d8d8d8', '{:}'.format(pdb1):'#000000'}

        df_plot = df_sorted.copy()
        f, ax = plt.subplots(1,1,figsize=(9,9))

        if df_dist.empty == False:
            ax.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])    
        
        df_plot = df_plot.drop(df_plot[df_plot.group == -1].index)
        df_other = df_plot.loc[df_plot["sort"] == 'noise']
        df_temp = df_other.rename({'i':'j', 'j':'i'}, axis=1)
        df_other = df_other.append(df_temp)
        df_plot = df_plot.drop(df_plot[df_plot.sort == 'noise'].index)
        df_c = df_plot.loc[df_plot["sort"] == 'common']
        df_c = df_c.rename({'i':'j', 'j':'i'}, axis=1)
        df_s = df_plot.loc[df_plot["sort"] == 'pdb_1']
        df_s = df_s.rename({'i':'j', 'j':'i'}, axis=1)
        df_plot = df_plot.drop(df_plot[df_plot.sort == 'pdb_1'].index)
        df_plot = df_plot.append(df_s)
        df_plot = df_plot.append(df_c)

        plot = ax.scatter(x=df_plot['i'], y=df_plot['j'], s=30, linewidth=0, c='#008080')
        
        df_other = df_other.drop_duplicates()
        ax.scatter(x=df_other['i'], y=df_other['j'],c='#008080', s=20,edgecolor='None',marker='D', alpha=0.3)

        #_______________________________________________________________________________________________________________________________
        #Final adjustments to image
        ax.set_xlabel('Residue Position', fontsize=30)
        ax.set_ylabel('Residue Position', fontsize=30)
        plt.xticks(fontsize=11, weight = 'bold')
        plt.yticks(fontsize=11, weight = 'bold')

        #_______________________________________________________________________________________________________________________________
        name = msa.split('/')
        plt.savefig('{:}.png'.format(name[-1][:-4]),dpi=600)
