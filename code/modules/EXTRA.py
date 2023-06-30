#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:42:44 2022

@author: schaferjw

Generate extra figures showing the difference in the top 1.5L contacts for the original MSA and smallest subfamily alignment 
of GREMLIN and MSATransformer. Generate figures showing all predictions, the effect of density based filtering, and bar graphs
showing the number of contacts corresponding to crstallographic strcutures.

"""
import pandas as pd
import numpy as np
from pathlib import Path
from os import listdir
from os.path import isfile, join
from modules.DB_SORT import DB_SCAN
from modules.SUPER import Data
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec
class GMN_Extra():
    def GMN_MSATR(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common): 
        #generate nondegenerate combination of gremlin and msa transformer predictions for comparison 
        gmn_file = Path('gmn/full_msa_00.csv')
        msatr_file = Path('msatr/full_msa_00.csv')

        if gmn_file.is_file() and msatr_file.is_file():
            df = pd.concat([self.df_original,self.df_original_msatr])
            df = df.drop_duplicates(subset=['i','j'])
            df['freq'] = 1
            
            self.df_gmn_msatr = df.copy()
            print('-------Results for Original GMN + MSATR-------')
            print('Total predictions: {:}'.format(len(df.index)))
            print('Percent within +-1 of first fold:  {:}'.format(len(df[df['sort'] == 'pdb_1'])))
            print('Percent within +-1 of second fold: {:}'.format(len(df[df['sort'] == 'pdb_2'])))
            print('Percent within +-1 of common fold: {:}'.format(len(df[df['sort'] == 'common'])))
            print('Average radial distance of noise from xcrystal : {:}'.format(df[df['sort'] == 'noise']['r'].mean()))
            print('\n')
      
    def DF_GMN_MSATR(self):return self.df_gmn_msatr
    def Super_No_Filter(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common,query_length):
	# generate the results of the superposition without filter
        df_super = pd.read_csv('{:}/df_super.csv'.format(msa[:-4]),index_col=0)
        if Prediction_Offset != 'n':
            df_super['i'] = df_super['i'] + int(Prediction_Offset)
            df_super['j'] = df_super['j'] + int(Prediction_Offset)

        db = DB_SCAN(df_super,xcontact,alignment)
        df_sorted = db.Run(x_first,x_second,x_common,msa,query_length, 'y')
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
    def MSA_Original_gmn(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common,query_length):
	# generate results for gremlin on the original msa
        gmn_file = Path('gmn/full_msa_00.csv')
        if gmn_file.is_file() :
            
            predictions = Data([gmn_file])
            df_super = predictions.Most_Probable(N_r=2.)
            
            if Prediction_Offset != 'n':
                df_super['i'] = df_super['i'] + int(Prediction_Offset)
                df_super['j'] = df_super['j'] + int(Prediction_Offset)
            df_super['freq'] = 1

            db = DB_SCAN(df_super,xcontact,alignment)
            df_sorted = db.Run(x_first,x_second,x_common,msa,query_length, 'y')
            self.df_original = df_sorted.copy()

            print('-------Results for Original GMN-------')
            print('Total predictions: {:}'.format(len(df_sorted.index)))
            print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
            print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
            print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
            print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
            print('\n')
    def MSA_Original_msatr(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common,query_length):
	#generate results for msa transformer on the original msa
        msatr_file = Path('msatr/full_msa_00.csv')
        if msatr_file.is_file():
            
            predictions = Data([msatr_file])
            df_super = predictions.Most_Probable(N_r=2.)
            
            if Prediction_Offset != 'n':
                df_super['i'] = df_super['i'] + int(Prediction_Offset)
                df_super['j'] = df_super['j'] + int(Prediction_Offset)
            df_super['freq'] = 1

            db = DB_SCAN(df_super,xcontact,alignment)
            df_sorted = db.Run(x_first,x_second,x_common,msa,query_length, 'y')
            self.df_original_msatr = df_sorted.copy()

            print('-------Results for Original MSATR-------')
            print('Total predictions: {:}'.format(len(df_sorted.index)))
            print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
            print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
            print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
            print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
            print('\n')
        
    def Subfamily_MSA_GMN(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common,query_length):
	# generate results of gremlin applied to the smallest subfamily alignment
        files = [i for i in listdir('gmn') if isfile(join('gmn', i))]
        files = [file for file in files if file[-4:] == '.csv']
        files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])     
        if files:
            gmn_file = Path('gmn/{:}'.format(files[-1]))
            if gmn_file.is_file():
                
                predictions = Data([gmn_file])
                df_super = predictions.Most_Probable(N_r=1.5)               
                
                if Prediction_Offset != 'n':
                    df_super['i'] = df_super['i'] + int(Prediction_Offset)
                    df_super['j'] = df_super['j'] + int(Prediction_Offset)
                self.sub = files[-1].split('_')[-1][:-4]
                df_super['freq'] = 1
                
                db = DB_SCAN(df_super,xcontact,alignment)
                df_sorted = db.Run(x_first,x_second,x_common,msa,query_length, 'y')
                self.df_sub = df_sorted.copy()

                print('-------Results for GMN sub-family: {:}-------'.format(files[-1]))
                print('Total predictions: {:}'.format(len(df_sorted.index)))
                print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
                print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
                print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
                print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
                print('\n')
    def Subfamily_MSA_MSATR(self,msa,xcontact, total_dist,Prediction_Offset,alignment,x_first,x_second,x_common,query_length):
	# generate results of msa transformer applied to the smallest subfamily alignment
        files = [i for i in listdir('msatr') if isfile(join('msatr', i))]
        files = [file for file in files if file[-4:] == '.csv']
        files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])     

        msatr_file = Path('msatr/{:}'.format(files[-1]))
        if msatr_file.is_file():

            predictions = Data([msatr_file])
            df_super = predictions.Most_Probable(N_r=1.5)               
 
            if Prediction_Offset != 'n':
                df_super['i'] = df_super['i'] + int(Prediction_Offset)
                df_super['j'] = df_super['j'] + int(Prediction_Offset)
            self.sub_msatr = files[-1].split('_')[-1][:-4]
            df_super['freq'] = 1
            
            db = DB_SCAN(df_super,xcontact,alignment)
            df_sorted = db.Run(x_first,x_second,x_common,msa,query_length, 'y')
            self.df_sub_msatr = df_sorted.copy()

            print('-------Results for MSATR sub-family: {:}-------'.format(files[-1]))
            print('Total predictions: {:}'.format(len(df_sorted.index)))
            print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_1'])))
            print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'pdb_2'])))
            print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[df_sorted['sort'] == 'common'])))
            print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted[df_sorted['sort'] == 'noise']['r'].mean()))
            print('\n')
    
    def SCATTER(self,ax,df,df_dist,colors,df_plot,plddt):
        if df_dist.empty == False:
            ax.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
        ax.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])
        
        
        if plddt[0] != 'none':
            alphas = ((0-1)/(np.max(plddt)-np.min(plddt)))*(plddt-np.min(plddt)) + 1
            alphas = [0.5 if alpha > 0.4 else 0.25 if alpha > 0.2 else 0 for alpha in alphas]

        else: alphas = plddt
        # scale for alphas is reversed so plddt of 80 - 100 gets no added color
        #                                 plddt of 60 - 80  gets yellow with alpha=0.25
        #                                 plddt of 0  - 60 gets red    with alpha=0.5
        plddt_rgba = [[143/255,20/255,2/255,alpha] if alpha == 0.5 else [255/255,191/255,0/255,alpha] if alpha == 0.25 else [0,0,0,0] for alpha in alphas]
        ax.scatter(x=df['i'], y=df['j'],c=plddt_rgba, s=100, linewidth=0)
        
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
        ax.set_xlabel('Residue Position', fontsize=20)
        ax.set_ylabel('Residue Position', fontsize=20)
        # ax.set_title('Original GMN', fontsize=25)
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=12)
        return plot
    
    def PLOT_EXTRA(self, msa, pdb1,pdb2,df_dist,df,df_sorted,plddt=['none']):

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
        if 'pdb_1'  not in counts: counts['pdb_1']  = 0
        if 'pdb_2'  not in counts: counts['pdb_2']  = 0
        colors = {'both':'#767676'}
        if counts['pdb_1'] > counts['pdb_2']:
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
        if counts['pdb_1'] <= counts['pdb_2']:
            colors = {'both':'#767676', '{:}'.format(pdb2):'#d8d8d8', '{:}'.format(pdb1):'#000000'}

         
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the original MSA
        #_______________________________________________________________________________________________________________________________
        plot = self.SCATTER(ax1,df,df_dist,colors,self.df_original,plddt)
        ax1.set_title('Original GMN', fontsize=25)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the shallowest sub-family alignment GMN
        #_______________________________________________________________________________________________________________________________
        ax3=fig.add_subplot(gs[1,0]) # Second row, third column
        ax_space2 = fig.add_subplot(gs[1,1])
        ax_space2.remove() #create space
        plot = self.SCATTER(ax3,df,df_dist,colors,self.df_sub,plddt)
        try:
            ax3.set_title('GMN Sub-Family: {:}'.format(self.sub), fontsize=25)
        except:
            self.sub = 'NA'
            ax3.set_title('GMN Sub-Family: {:}'.format(self.sub), fontsize=25)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of MSATR original
        #_______________________________________________________________________________________________________________________________
        ax4=fig.add_subplot(gs[0,2]) # Second row, First column
        ax_space3 = fig.add_subplot(gs[0,3])
        ax_space3.remove() #create space
        plot = self.SCATTER(ax4,df,df_dist,colors,self.df_original_msatr,plddt)
        ax4.set_title('Original MSATR', fontsize=25)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the shallowest sub-family alignment MSATR
        #_______________________________________________________________________________________________________________________________
        ax5=fig.add_subplot(gs[1,2]) # Second row, third column
        ax_space4 = fig.add_subplot(gs[1,3])
        ax_space4.remove() #create space
        plot = self.SCATTER(ax5,df,df_dist,colors,self.df_sub_msatr,plddt)
        ax5.set_title('MSATR Sub-Family: {:}'.format(self.sub_msatr), fontsize=25)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the full pipeline
        #_______________________________________________________________________________________________________________________________
        ax7=fig.add_subplot(gs[0,4]) # First row, first column
        plot = self.SCATTER(ax7,df,df_dist,colors,df_sorted.drop(df_sorted[df_sorted.group == -1].index),plddt)
        df_plot = df_sorted.drop(df_sorted[df_sorted.group == -1].index)
        df_full = df_plot.copy()
        ax7.set_title('Full Pipeline', fontsize=25)
        
        #_______________________________________________________________________________________________________________________________
        # Plot the result of the unfiltered superposition
        #_______________________________________________________________________________________________________________________________
        ax9=fig.add_subplot(gs[1,4]) # First row, third column
        plot = self.SCATTER(ax9,df,df_dist,colors,self.df_no_filter,plddt)
        ax9.set_title('Superposition', fontsize=25)

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
            
            if counts['pdb_1'] <= counts['pdb_2']:
                rects1 = ax11.bar(df_stack.index, df_stack.Noise,  width, color='#008080',label='Noise')
                rects2 = ax11.bar(df_stack.index, df_stack.Common, width, bottom=df_stack.Noise, color='#767676',label='common')
                rects3 = ax11.bar(df_stack.index, df_stack.First,  width, bottom=df_stack.Noise+df_stack.Common, color='#000000',label=pdb1)
                rects4 = ax11.bar(df_stack.index, df_stack.Second, width, bottom=df_stack.Noise+df_stack.Common+df_stack.First, color='#d8d8d8',label=pdb2)
            
            if counts['pdb_1'] > counts['pdb_2']:
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
        plt.savefig('Extra_{:}.png'.format(name[-1][:-4]),dpi=200)
