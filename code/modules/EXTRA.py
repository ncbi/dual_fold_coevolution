#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:42:44 2022

@author: schaferjw
"""
import pandas as pd
import numpy as np
import scipy
from pathlib import Path
from os import listdir
from os.path import isfile, join
from modules.DB_SORT import DB_SCAN
import matplotlib
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec
class GMN_Extra():
    def GMN_MSATR(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common,N): 
        #generate nondegenerate combination of gremlin and msa transformer predictions for comparison 
        gmn_file = Path('{:}/df_gmn_00.csv'.format(msa[:-4]))
        msatr_file = Path('{:}/df_msatr_00.csv'.format(msa[:-4]))
        if gmn_file.is_file() and msatr_file.is_file():
            df_gmn = pd.read_csv('{:}/df_gmn_00.csv'.format(msa[:-4]),index_col=0)
            df_msatr = pd.read_csv('{:}/df_msatr_00.csv'.format(msa[:-4]),index_col=0)
            df = pd.concat([df_gmn,df_msatr])
            df = df.drop_duplicates(subset=['i','j'])
            df['freq'] = 1
            if Prediction_Offset != 'n':
                df['i'] = df['i'] + int(Prediction_Offset)
                df['j'] = df['j'] + int(Prediction_Offset)
            db = DB_SCAN(df,xcontact,alignment)
            df_sorted = db.Run(x_first,x_second,x_common,msa, 'y')

            self.df_gmn_msatr = df_sorted        
            print('-------Results for Original GMN + MSATR-------')
            print('Total predictions: {:}'.format(len(df_sorted.index)))
            print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
            print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
            print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
            print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
            print('\n')
        
    def DF_GMN_MSATR(self):return self.df_gmn_msatr
    def Super_No_Filter(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common,N):
	# generate the results of the superposition without filter
        df_super = pd.read_csv('{:}/df_super.csv'.format(msa[:-4]),index_col=0)
        if Prediction_Offset != 'n':
            df_super['i'] = df_super['i'] + int(Prediction_Offset)
            df_super['j'] = df_super['j'] + int(Prediction_Offset)

        db = DB_SCAN(df_super,xcontact,alignment)
        df_sorted = db.Run(x_first,x_second,x_common,msa, 'y')
        self.df_no_filter = df_sorted.copy()
        
        print('-------Results for full Pipeline no filter-------')
        print('Total predictions: {:}'.format(len(df_sorted.index)))
        print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
        print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
        print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
        print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
        print('\n')

    def DF_no_filter(self): return self.df_no_filter
    def DF_original(self): return self.df_original
    def DF_orig_msatr(self): return self.df_original_msatr
    def MSA_Original_gmn(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common):
	# generate results for gremlin on the original msa
        gmn_file = Path('{:}/df_gmn_00.csv'.format(msa[:-4]))
        if gmn_file.is_file() :
            df_super = pd.read_csv('{:}/df_gmn_00.csv'.format(msa[:-4]),index_col=0)
            if Prediction_Offset != 'n':
                df_super['i'] = df_super['i'] + int(Prediction_Offset)
                df_super['j'] = df_super['j'] + int(Prediction_Offset)
            df_super['freq'] = 1

            db = DB_SCAN(df_super,xcontact,alignment)
            df_sorted = db.Run(x_first,x_second,x_common,msa, 'y')
            self.df_original = df_sorted.copy()

            print('-------Results for Original GMN-------')
            print('Total predictions: {:}'.format(len(df_sorted.index)))
            print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
            print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
            print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
            print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
            print('\n')
    def MSA_Original_msatr(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common):
	#generate results for msa transformer on the original msa
        msatr_file = Path('{:}/df_msatr_00.csv'.format(msa[:-4]))
        if msatr_file.is_file():
            df_super = pd.read_csv('{:}/df_msatr_00.csv'.format(msa[:-4]),index_col=0)
            if Prediction_Offset != 'n':
                df_super['i'] = df_super['i'] + int(Prediction_Offset)
                df_super['j'] = df_super['j'] + int(Prediction_Offset)
            df_super['freq'] = 1

            db = DB_SCAN(df_super,xcontact,alignment)
            df_sorted = db.Run(x_first,x_second,x_common,msa, 'y')
            self.df_original_msatr = df_sorted.copy()

            print('-------Results for Original MSATR-------')
            print('Total predictions: {:}'.format(len(df_sorted.index)))
            print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
            print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
            print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
            print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
            print('\n')
        
    def Subfamily_MSA_GMN(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common):
	# generate results of gremlin applied to the smallest subfamily alignment
        files = [i for i in listdir('{:}'.format(msa[:-4])) if isfile(join('{:}'.format(msa[:-4]), i))]
        temp = [file for file in files if file[:6] == 'df_gmn' and file != 'df_super.csv']
        files = temp
        files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])     
        if files:
            gmn_file = Path('{:}/{:}'.format(msa[:-4],files[-1]))
            if gmn_file.is_file():

                df_super = pd.read_csv('{:}/{:}'.format(msa[:-4],files[-1]),index_col=0)
                if Prediction_Offset != 'n':
                    df_super['i'] = df_super['i'] + int(Prediction_Offset)
                    df_super['j'] = df_super['j'] + int(Prediction_Offset)
                self.sub = files[-1].split('_')[2][:-4]
                df_super['freq'] = 1
                
                db = DB_SCAN(df_super,xcontact,alignment)
                df_sorted = db.Run(x_first,x_second,x_common,msa, 'y')
                self.df_sub = df_sorted.copy()

                print('-------Results for GMN sub-family: {:}-------'.format(files[-1]))
                print('Total predictions: {:}'.format(len(df_sorted.index)))
                print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
                print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
                print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
                print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
                print('\n')
    def Subfamily_MSA_MSATR(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common):
	# generate results of msa transformer applied to the smallest subfamily alignment
        files = [i for i in listdir('{:}'.format(msa[:-4])) if isfile(join('{:}'.format(msa[:-4]), i))]
        files = [file for file in files if file[:8] == 'df_msatr' and file != 'df_super.csv']
        files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])     

        msatr_file = Path('{:}/{:}'.format(msa[:-4],files[-1]))
        if msatr_file.is_file():

            df_super = pd.read_csv('{:}/{:}'.format(msa[:-4],files[-1]),index_col=0)
            if Prediction_Offset != 'n':
                df_super['i'] = df_super['i'] + int(Prediction_Offset)
                df_super['j'] = df_super['j'] + int(Prediction_Offset)
            self.sub_msatr = files[-1].split('_')[2][:-4]
            df_super['freq'] = 1
            
            db = DB_SCAN(df_super,xcontact,alignment)
            df_sorted = db.Run(x_first,x_second,x_common,msa, 'y')
            self.df_sub_msatr = df_sorted.copy()

            print('-------Results for MSATR sub-family: {:}-------'.format(files[-1]))
            print('Total predictions: {:}'.format(len(df_sorted.index)))
            print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
            print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
            print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
            print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
            print('\n')
    
    def PLOT_EXTRA(self, msa, pdb1,pdb2,df_dist,df,df_sorted):

        cmap = plt.get_cmap("plasma")
        fig=plt.figure(figsize=(24,20))
        
        gs=GridSpec(3,5, width_ratios=[1,0.35,1,0.35,1], height_ratios=[1,1,1.]) # 5 rows, 7 columns
        gs.update(wspace=0., hspace=0.25)
        
        #_______________________________________________________________________________________________________________________________
        # load all alignments needed for generating the 6 panel png with bar graph
        #_______________________________________________________________________________________________________________________________
        ax1=fig.add_subplot(gs[0,0]) # Second row, First column
        ax_space = fig.add_subplot(gs[0,1])
        ax_space.remove() #create space
        
        counts = df_sorted['sort'].value_counts()
        colors = {'both':'#767676'}
        if len(counts) == 4 and counts['pdb_1'] > counts['pdb_2']:
            df_sorted = df_sorted.rename({'i':'j', 'j':'i'}, axis=1)
            df_dist = df_dist.rename({'i':'j', 'j':'i'}, axis=1)
            df = df.rename({'i':'j', 'j':'i'}, axis=1)
            try:
                self.df_original = self.df_original.rename({'i':'j', 'j':'i'}, axis=1)
            except:
                self.df_original = pd.DataFrame(columns=["i","j","apc","zscore","gmn_x_freq","group","sort","r","df","freq"])
            try:
                self.df_sub = self.df_sub.rename({'i':'j', 'j':'i'}, axis=1)
            except:
                self.df_sub = pd.DataFrame(columns=["i","j","apc","zscore","gmn_x_freq","group","sort","r","df","freq"])
            try:
                self.df_original_msatr = self.df_original_msatr.rename({'i':'j', 'j':'i'}, axis=1)
            except:
                self.df_original_msatr = pd.DataFrame(columns=["i","j","apc","zscore","gmn_x_freq","group","sort","r","df","freq"])
            try:
                self.df_sub_msatr = self.df_sub_msatr.rename({'i':'j', 'j':'i'}, axis=1)
            except:
                self.df_sub_msatr = pd.DataFrame(columns=["i","j","apc","zscore","gmn_x_freq","group","sort","r","df","freq"])
            self.df_no_filter = self.df_no_filter.rename({'i':'j', 'j':'i'}, axis=1)
            #Make sure both is always gray
            colors = {'both':'#767676', '{:}'.format(pdb1):'#d8d8d8', '{:}'.format(pdb2):'#000000'}
        if len(counts) == 4 and counts['pdb_1'] < counts['pdb_2']:
            colors = {'both':'#767676', '{:}'.format(pdb2):'#d8d8d8', '{:}'.format(pdb1):'#000000'}

         
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the original MSA
        #_______________________________________________________________________________________________________________________________
    
        if df_dist.empty == False:
            ax1.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax1.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])
        
        df_plot = self.df_original
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
  
        plot = ax1.scatter(x=df_plot['i'], y=df_plot['j'], s=30, linewidth=0, c='#008080')
        df_other = df_other.drop_duplicates()
        ax1.scatter(x=df_other['i'], y=df_other['j'],c='#008080', s=20,edgecolor='None',marker='D', alpha=0.3)
        ax1.set_xlabel('Residue Position', fontsize=20)
        ax1.set_ylabel('Residue Position', fontsize=20)
        ax1.set_title('Original GMN', fontsize=25)
        ax1.tick_params(axis="x", labelsize=12)
        ax1.tick_params(axis="y", labelsize=12)
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the shallowest sub-family alignment GMN
        #_______________________________________________________________________________________________________________________________
        ax3=fig.add_subplot(gs[1,0]) # Second row, third column
        ax_space2 = fig.add_subplot(gs[1,1])
        ax_space2.remove() #create space
        
        if df_dist.empty == False:
            ax3.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax3.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])
        
        df_plot = self.df_sub
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
       
        plot = ax3.scatter(x=df_plot['i'], y=df_plot['j'], s=30, linewidth=0, c='#008080')
        
        df_other = df_other.drop_duplicates()
       
        ax3.scatter(x=df_other['i'], y=df_other['j'],c='#008080', s=20,edgecolor='None',marker='D', alpha=0.3)
        ax3.set_xlabel('Residue Position', fontsize=20)
        ax3.set_ylabel('Residue Position', fontsize=20)
        try:
            ax3.set_title('GMN Sub-Family: {:}'.format(self.sub), fontsize=25)
        except:
            self.sub = 'NA'
            ax3.set_title('GMN Sub-Family: {:}'.format(self.sub), fontsize=25)
        ax3.tick_params(axis="x", labelsize=12)
        ax3.tick_params(axis="y", labelsize=12)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of MSATR original
        #_______________________________________________________________________________________________________________________________
        ax4=fig.add_subplot(gs[0,2]) # Second row, First column
        ax_space3 = fig.add_subplot(gs[0,3])
        ax_space3.remove() #create space
        if df_dist.empty == False:
            ax4.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax4.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])
        
        df_plot = self.df_original_msatr
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
  
        plot = ax4.scatter(x=df_plot['i'], y=df_plot['j'], s=30, linewidth=0, c='#008080')
        df_other = df_other.drop_duplicates()
  
        ax4.scatter(x=df_other['i'], y=df_other['j'],c='#008080', s=20,edgecolor='None',marker='D', alpha=0.3)
        ax4.set_xlabel('Residue Position', fontsize=20)
        ax4.set_ylabel('Residue Position', fontsize=20)
        ax4.set_title('Original MSATR', fontsize=25)
        ax4.tick_params(axis="x", labelsize=12)
        ax4.tick_params(axis="y", labelsize=12)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the shallowest sub-family alignment MSATR
        #_______________________________________________________________________________________________________________________________
        ax5=fig.add_subplot(gs[1,2]) # Second row, third column
        ax_space4 = fig.add_subplot(gs[1,3])
        ax_space4.remove() #create space
        
        if df_dist.empty == False:
            ax5.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax5.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])
        
        df_plot = self.df_sub_msatr
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
       
        plot = ax5.scatter(x=df_plot['i'], y=df_plot['j'], s=30, linewidth=0, c='#008080')
        df_other = df_other.drop_duplicates()
        ax5.scatter(x=df_other['i'], y=df_other['j'],c='#008080', s=20,edgecolor='None',marker='D', alpha=0.3)

        ax5.set_xlabel('Residue Position', fontsize=20)
        ax5.set_ylabel('Residue Position', fontsize=20)
        ax5.set_title('MSATR Sub-Family: {:}'.format(self.sub_msatr), fontsize=25)
        ax5.tick_params(axis="x", labelsize=12)
        ax5.tick_params(axis="y", labelsize=12)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the full pipeline
        #_______________________________________________________________________________________________________________________________

        ax7=fig.add_subplot(gs[0,4]) # First row, first column
    
        if df_dist.empty == False:
            ax7.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax7.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])    
        df_plot = df_sorted.drop(df_sorted[df_sorted.group == -1].index)
        df_full = df_plot.copy()
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

        plot = ax7.scatter(x=df_plot['i'], y=df_plot['j'], s=30, linewidth=0, c='#008080')
        df_other = df_other.drop_duplicates()
        ax7.scatter(x=df_other['i'], y=df_other['j'],c='#008080', s=20,edgecolor='None',marker='D', alpha=0.3)
       
        #_______________________________________________________________________________________________________________________________
        #Final adjustments to image
        ax7.set_xlabel('Residue Position', fontsize=20)
        ax7.set_ylabel('Residue Position', fontsize=20)
        ax7.set_title('Full Pipeline', fontsize=25)
        ax7.tick_params(axis="x", labelsize=12)
        ax7.tick_params(axis="y", labelsize=12)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the unfiltered superposition
        #_______________________________________________________________________________________________________________________________
        ax9=fig.add_subplot(gs[1,4]) # First row, third column
        
        if df_dist.empty == False:
            ax9.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax9.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])
        
        df_plot = self.df_no_filter
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
  
        plot = ax9.scatter(x=df_plot['i'], y=df_plot['j'], s=30, linewidth=0, c='#008080')
        df_other = df_other.drop_duplicates()
        ax9.scatter(x=df_other['i'], y=df_other['j'],c='#008080', s=20,edgecolor='None',marker='D', alpha=0.3)
       

        ax9.set_xlabel('Residue Position', fontsize=20)
        ax9.set_ylabel('Residue Position', fontsize=20)
        ax9.set_title('Superposition', fontsize=25)
        ax9.tick_params(axis="x", labelsize=12)
        ax9.tick_params(axis="y", labelsize=12)
      
        #_______________________________________________________________________________________________________________________________
        # Bar graph of changes in distribution of predictions
        #_______________________________________________________________________________________________________________________________
        
        if len(self.df_sub) != 0:
            ax11=fig.add_subplot(gs[2,:]) # Third row, All columns
            
            data_dict = {'Original GMN' : [len(self.df_original[self.df_original['sort'] == 'noise'].index),
                        len(self.df_original[self.df_original['sort'] == 'common'].index),
                        len(self.df_original[self.df_original['sort'] == 'pdb_1'].index),
                        len(self.df_original[self.df_original['sort'] == 'pdb_2'].index)],
            
            'Sub-Family GMN: {:}'.format(self.sub): [len(self.df_sub[self.df_sub['sort'] == 'noise'].index),
                        len(self.df_sub[self.df_sub['sort'] == 'common'].index),
                        len(self.df_sub[self.df_sub['sort'] == 'pdb_1'].index),
                        len(self.df_sub[self.df_sub['sort'] == 'pdb_2'].index)],
            
            'Original MSATR' : [len(self.df_original_msatr[self.df_original_msatr['sort'] == 'noise'].index),
                        len(self.df_original_msatr[self.df_original_msatr['sort'] == 'common'].index),
                        len(self.df_original_msatr[self.df_original_msatr['sort'] == 'pdb_1'].index),
                        len(self.df_original_msatr[self.df_original_msatr['sort'] == 'pdb_2'].index)],
            
            'Sub-Family MSATR: {:}'.format(self.sub_msatr): [len(self.df_sub_msatr[self.df_sub_msatr['sort'] == 'noise'].index),
                        len(self.df_sub_msatr[self.df_sub_msatr['sort'] == 'common'].index),
                        len(self.df_sub_msatr[self.df_sub_msatr['sort'] == 'pdb_1'].index),
                        len(self.df_sub_msatr[self.df_sub_msatr['sort'] == 'pdb_2'].index)],
            
            'Full Pipeline': [len(df_full[df_full['sort'] == 'noise'].index),
                        len(df_full[df_full['sort'] == 'common'].index),
                        len(df_full[df_full['sort'] == 'pdb_1'].index),
                        len(df_full[df_full['sort'] == 'pdb_2'].index)],
            
            'Superposition': [len(self.df_no_filter[self.df_no_filter['sort'] == 'noise'].index),
                        len(self.df_no_filter[self.df_no_filter['sort'] == 'common'].index),
                        len(self.df_no_filter[self.df_no_filter['sort'] == 'pdb_1'].index),
                        len(self.df_no_filter[self.df_no_filter['sort'] == 'pdb_2'].index)]}
            
            labels = ['Noise', 'Common', 'First', 'Second']
            cmap = plt.get_cmap("inferno")
            
            df_stack = pd.DataFrame.from_dict(data_dict, orient='index', columns=labels)        
            ax11.set_ylabel('Number of Predictions',fontsize=20)
            ax11.tick_params(axis="x", labelsize=20)
            ax11.tick_params(axis="y", labelsize=20)

            width = 0.35
            
            #default for single fold protiens
            if len(counts) == 2:
                rects1 = ax11.bar(df_stack.index, df_stack.Noise,  width, color='#008080',label='Noise')
                rects2 = ax11.bar(df_stack.index, df_stack.Common, width, bottom=df_stack.Noise, color='#767676',label='common')
                rects3 = ax11.bar(df_stack.index, df_stack.First,  width, bottom=df_stack.Noise+df_stack.Common, color='#d8d8d8',label=pdb1)
                rects4 = ax11.bar(df_stack.index, df_stack.Second, width, bottom=df_stack.Noise+df_stack.Common+df_stack.First, color='#000000',label=pdb2)
            
            if len(counts) == 4 and counts['pdb_1'] < counts['pdb_2']:
                rects1 = ax11.bar(df_stack.index, df_stack.Noise,  width, color='#008080',label='Noise')
                rects2 = ax11.bar(df_stack.index, df_stack.Common, width, bottom=df_stack.Noise, color='#767676',label='common')
                rects3 = ax11.bar(df_stack.index, df_stack.First,  width, bottom=df_stack.Noise+df_stack.Common, color='#000000',label=pdb1)
                rects4 = ax11.bar(df_stack.index, df_stack.Second, width, bottom=df_stack.Noise+df_stack.Common+df_stack.First, color='#d8d8d8',label=pdb2)
            
            if len(counts) == 4 and counts['pdb_1'] > counts['pdb_2']:
                rects1 = ax11.bar(df_stack.index, df_stack.Noise,  width, color='#008080',label='Noise')
                rects2 = ax11.bar(df_stack.index, df_stack.Common, width, bottom=df_stack.Noise, color='#767676',label='common')
                rects3 = ax11.bar(df_stack.index, df_stack.First,  width, bottom=df_stack.Noise+df_stack.Common, color='#d8d8d8',label=pdb1)
                rects4 = ax11.bar(df_stack.index, df_stack.Second, width, bottom=df_stack.Noise+df_stack.Common+df_stack.First, color='#000000',label=pdb2)
            
            ax11.legend()

            for r1, r2, r3, r4 in zip(rects1, rects2, rects3, rects4):
                h1 = r1.get_height()
                h2 = r2.get_height()
                h3 = r3.get_height()
                h4 = r4.get_height()
                
                percent = np.around(int(h1 * 100 / (h1+h2+h3+h4)),decimals=0)
                ax11.annotate('{}%'.format(percent),   #percent for noise
                                xy=(r1.get_x() - r1.get_width() / 5, h1/2),
                                xytext=(0, 0),
                                textcoords="offset points",
                                ha='center', va='bottom', color='k', weight='bold')
                ax11.annotate('{}'.format(r1.get_height()),  #number for noise
                                xy=(r1.get_x() + r1.get_width() +r1.get_width()/5, h1/2),
                                xytext=(0, 0),
                                textcoords="offset points",
                                ha='center', va='bottom', color='k', weight='bold')
                percent = np.around(int(h2 * 100 / (h1+h2+h3+h4)),decimals=0)
                ax11.annotate('{}%'.format(percent),   #percent for common
                                xy=(r2.get_x() - r2.get_width() / 5, h1+h2/2),
                                xytext=(0, 0), 
                                textcoords="offset points",
                                ha='center', va='bottom', color='k', weight='bold')
                ax11.annotate('{}'.format(r2.get_height()),  #number for common
                                xy=(r2.get_x() + r2.get_width() +r2.get_width()/5, h1+h2/2),
                                xytext=(0, 0),
                                textcoords="offset points",
                                ha='center', va='bottom', color='k', weight='bold')
                percent = np.around(int(h3 * 100 / (h1+h2+h3+h4)),decimals=0)
                ax11.annotate('{}%'.format(percent),   #percent for first
                                xy=(r3.get_x() - r3.get_width() / 5, h1+h2+h3/2),
                                xytext=(0, 0), 
                                textcoords="offset points",
                                ha='center', va='bottom', color='k', weight='bold')
                ax11.annotate('{}'.format(r3.get_height()),   #number for first
                                xy=(r3.get_x() + r3.get_width() +r3.get_width()/5, h1+h2+h3/2),
                                xytext=(0, 0),
                                textcoords="offset points",
                                ha='center', va='bottom', color='k', weight='bold')
                percent = np.around(int(h4 * 100 / (h1+h2+h3+h4)),decimals=0)
                ax11.annotate('{}%'.format(percent),   #percent for second
                                xy=(r4.get_x() - r4.get_width() / 5, h1+h2+h3+h4/2),
                                xytext=(0, 0), 
                                textcoords="offset points",
                                ha='center', va='bottom', color='k', weight='bold')
                ax11.annotate('{}'.format(r4.get_height()),   #number for second
                                xy=(r4.get_x() + r4.get_width() +r4.get_width()/5, h1+h2+h3+h4/2),
                                xytext=(0, 0),
                                textcoords="offset points",
                                ha='center', va='bottom', color='k', weight='bold')
            #_______________________________________________________________________________________________________________________________

        name = msa.split('/')
        plt.savefig('Extra_{:}.png'.format(name[-1][:-4]),dpi=400)
