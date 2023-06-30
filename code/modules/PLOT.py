#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 15:58:10 2022

@author: schaferjw

Plot predictions on top of the dualfold contact map

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class PLOT():
    def FS_plot(self,pdb1,pdb2,df_dist,df,msa,df_sorted,name='n',plddt=['none']):
        colors = {'both':'#767676'}
        counts = df_sorted['sort'].value_counts()
        df['plddt'] = plddt
        if 'pdb_1' not in counts: counts['pdb_1']  = 0
        if 'pdb_2' not in counts: counts['pdb_2']  = 0
        if counts['pdb_1'] > counts['pdb_2']:
            df_sorted = df_sorted.rename({'i':'j', 'j':'i'}, axis=1)
            df_dist = df_dist.rename({'i':'j', 'j':'i'}, axis=1)
            df = df.rename({'i':'j', 'j':'i'}, axis=1)
            colors = {'both':'#767676', '{:}'.format(pdb1):'#d8d8d8', '{:}'.format(pdb2):'#000000'}
        if counts['pdb_1'] <= counts['pdb_2']:
            colors = {'both':'#767676', '{:}'.format(pdb2):'#d8d8d8', '{:}'.format(pdb1):'#000000'}

        if plddt[0] != 'none':
            alphas = ((0-1)/(np.max(plddt)-np.min(plddt)))*(plddt-np.min(plddt)) + 1
            alphas = [0.5 if alpha > 0.4 else 0.25 if alpha > 0.2 else 0 for alpha in alphas]

        else: alphas = plddt
        # scale for alphas is reversed so plddt of 80 - 100 gets no added color
        #                                 plddt of 60 - 80  gets yellow with alpha=0.25
        #                                 plddt of 0  - 60  gets red    with alpha=0.5
        plddt_rgba = [[143/255,20/255,2/255,alpha] if alpha == 0.5 else [255/255,191/255,0/255,alpha] if alpha == 0.25 else [0,0,0,0] for alpha in alphas]
        df_plot = df_sorted.copy()
        f, ax = plt.subplots(1,1,figsize=(9,9))

        if df_dist.empty == False:
            ax.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])
        ax.scatter(x=df['i'], y=df['j'],c=plddt_rgba, s=100, linewidth=0)
        
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
        if name == 'n':
            name = msa.split('/')
            plt.savefig('{:}.png'.format(name[-1][:-4]),dpi=200)
        else:
            plt.savefig('{:}.png'.format(name),dpi=200)
            