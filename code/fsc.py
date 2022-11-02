#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This code is intended for evaluating the overlap between crystallographic contacts of known fold-switching proteins
and coevolutionary predictions made using GREMLIN and MSA Transformer.

To run this code .pdb files are required for each conformation and an MSA file is needed in stockholm format.
NOTE: the stockholm files were generated with HMMER 3.3.2 and renamed as 'pdbID.msa' (this program expects the .msa extension
@author: schaferjw
"""

import argparse
import pandas as pd
import numpy as np
import os
import sys
import filecmp
from os import listdir
from pathlib import Path
from os.path import isfile, join

from modules.PDB import PDB_STRUCTURE
from modules.MSA_PREP import HHFILTER
from modules.MSA_TR import MSA_TR
from modules.COEVOLUTION import Start
from modules.SUPER import GMN_Data
from modules.DB_SORT import DB_SCAN
from modules.PLOT import PLOT
from modules.EXTRA import GMN_Extra
from modules.HYPERGEO import HYPERGEO

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--msa", type=str, help='MSA file in stockholm format')
    parser.add_argument("--pdb1", type=str, help='PDB structure for the first conformation (remove hydrogens, ligands, and HOH')
    parser.add_argument("--pdb2", type=str, help='PDB structure for the second conformation (remove hydrogens, ligands, and HOH')
    parser.add_argument("--Extra", type=str,default='n', help='Output additional information on superposition, base case, and smallest subfamily alignment')
    parser.add_argument("--Depth_max", type=float,default=2.0, help='Control how many predictions are extracted from subfamily alignments (larger number = more predictions)')
    parser.add_argument("--Depth_min", type=float,default=1.5, help='Control how many predictions are extracted from subfamily alignments (larger number = more predictions)')
    parser.add_argument("--Manual_PDB", type=str,default='n', help='Turn off auto-align algorithm and manually align PDBs')
    parser.add_argument("--Manual_Pred", type=str,default='n', help='Turn off auto-align algorithm and manually align Predictions')
    parser.add_argument("--Beg", type=int,default=0, help='Remove x number of residues from the beginning of the msa')
    parser.add_argument("--End", type=int,default=0, help='Remove x number of residues from the end of the msa')
    parser.add_argument("--Prediction_Offset", type=str,default='n', help='Used if a large number of residues are missing from the CTD (ex. 3ejh/3m7p alignment set to 225')
    
    args = parser.parse_args()
    
    
    #   *** NOTE: This section will not run if the directory it makes is already present ***
    #   Prepare Subfamily MSAs from alignment generated from HMMER
    #___________________________________________________________________________________________________
    
    if os.path.isdir(args.msa[:-4]) == False:
        prep = HHFILTER(args.msa)
        prep.Filter()
        prep.Data(args.Beg,args.End)
    #___________________________________________________________________________________________________
    
    #   Run msa transformer
    #___________________________________________________________________________________________________
    
    
        files = [i for i in listdir('{:}/'.format(args.msa[:-4])) if isfile(join('{:}/'.format(args.msa[:-4]), i))]
        files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])
        remove = []
        files = [file for file in files if file[-4:] == '.a3m']
        for i in files:
            num_temp = sum(1 for line in open('{:}/{:}'.format(args.msa[:-4],i)))
            if num_temp/2 < 8 and i[0:2] != '00': 
                os.system('rm -f {:}/{:}'.format(args.msa[:-4],i))
                remove.append(i)
        files_tr = [file for file in files if file not in remove]
        remove = []
        for i in range(len(files_tr)-1):
            if filecmp.cmp('{:}/{:}'.format(args.msa[:-4],files_tr[i]), '{:}/{:}'.format(args.msa[:-4],files_tr[i+1])) == True:
                remove.append(files_tr[i+1])
        for i in remove:
            os.system('rm -f {:}/{:}'.format(args.msa[:-4],i))
        files_tr = [file for file in files_tr if file not in remove]
        MSA_TR(args.msa[:-4],files_tr,args.Beg,args.End)
    
    #___________________________________________________________________________________________________
    #   Run Coevolutionary analysis
    #___________________________________________________________________________________________________
    
        files = [i for i in listdir('{:}/'.format(args.msa[:-4])) if isfile(join('{:}/'.format(args.msa[:-4]), i))]
        files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])
        files_gmn = [file for file in files if file[-4:] == '.out']
        with open('{:}/{:}'.format(args.msa[:-4],files_gmn[0])) as f:
            lines = [line.rstrip() for line in f]
            cutoff = len(lines[0].split(' ')[-1:][0])
        
        print('L : {:}'.format(cutoff))
        print('msa depth: {:}'.format(len(lines))) 
    
        remove = []
        for i in files_gmn:
            num_temp = sum(1 for line in open('{:}/{:}'.format(args.msa[:-4],i)))
            if num_temp < 5*cutoff: 
                os.system('rm -f {:}/{:}'.format(args.msa[:-4],i))
                os.system('rm -f {:}/{:}.csv'.format(args.msa[:-4],i[:-4]))
                os.system('rm -f {:}/{:}.sto'.format(args.msa[:-4],i[:-4]))
                remove.append(i)
        files_gmn = [file for file in files_gmn if file not in remove]
        remove = []
        for i in range(len(files_gmn)-1): #remove redundant msas
            if filecmp.cmp('{:}/{:}'.format(args.msa[:-4],files_gmn[i]), '{:}/{:}'.format(args.msa[:-4],files_gmn[i+1])) == True:
                remove.append(files_gmn[i+1])
        files_gmn = [file for file in files_gmn if file not in remove]
        for i in remove:
            os.system('rm -f {:}/{:}'.format(args.msa[:-4],i))
            os.system('rm -f {:}/{:}.csv'.format(args.msa[:-4],i[:-4]))
            os.system('rm -f {:}/{:}.sto'.format(args.msa[:-4],i[:-4]))
    
        output = np.array([])
        if not files_gmn:
            print('...Not enough data for GREMLIN algorithm...')
        else:
            for i in files_gmn: #scale prediction depth based on change in msa depth
                num_temp = sum(1 for line in open('{:}/{:}'.format(args.msa[:-4],i)))
                output = np.append(output,num_temp)
            output = ((args.Depth_max-args.Depth_min)/(np.max(output)-np.min(output)))*(output-np.min(output)) + args.Depth_min
            output[::-1].sort()
            if np.isnan(output[0]) == True:
                output[0] = 2.0
            for file in zip(files_gmn,np.reshape(output,(len(files_gmn),1)).tolist()):
                print('Coevolutionary Analysis running on {:}'.format(file[0]))
                prd = Start(args.msa[:-4],file[0],file[1][0])
                prd.to_csv('{:}/df_gmn_{:}.csv'.format(args.msa[:-4],file[0][:-4]))
    
    
    #___________________________________________________________________________________________________
    #   Superimpose predictions into one contact map
    #___________________________________________________________________________________________________
    
        predictions = GMN_Data(args.msa)
        predictions.Super()
    
    
    #___________________________________________________________________________________________________
    
    #   Load crystallagraphic information
    #___________________________________________________________________________________________________
    #MDtraj uses deprecated np.int instead of the new np.int32 or np.int64
    #comment these lines out to see full deprecation warning
    import warnings
    warnings.filterwarnings('ignore')
    
    pdbs = PDB_STRUCTURE(args.pdb1,args.pdb2)
    pdbs.Get_contacts(1.)       #define list of intrachain contacts for comparison
    pdbs.Multimer_contact(1.)   #create list of inerchain contacts (based on above)
    pdbs.Get_contacts(0.8)      #redefine list of intrachain contacts at 8A
    pdbs.Combo_opt()
    fold_id, dist_fold_id, opt, pdb_1, pdb_2, dist_pdb_1, dist_pdb_2 = pdbs.Combo(args.Manual_PDB,args.msa[:-4])
    xcontact = pdb_1 + pdb_2
    df = pd.DataFrame(xcontact,columns=list('ij'))
    df['Fold'] = fold_id
    total_dist = dist_pdb_1 + dist_pdb_2
    df_dist = pd.DataFrame(total_dist,columns=list('ij'))
    df_dist['Fold'] = dist_fold_id
    
    #___________________________________________________________________________________________________
    #   Find Clusters within predictions
    #   align superposition of predictions to crystallographic information and cluster
    #   Print information about xcrystal alignment, predictions, and save png file.
    #___________________________________________________________________________________________________
    if os.path.isdir(args.msa[:-4]) == True:print('Precomputed Coevolutionary predictions are being loaded from driectory {:}/ '.format(args.msa[:-4]))
    df_super = pd.read_csv('{:}/df_super.csv'.format(args.msa[:-4]),index_col=0)
    if args.Prediction_Offset != 'n':
        df_super['i'] = df_super['i'] + int(args.Prediction_Offset)
        df_super['j'] = df_super['j'] + int(args.Prediction_Offset)
    
    db = DB_SCAN(df_super,xcontact,args.Manual_Pred)
    alignment = db.Return_Alignment()
    
    df_tot = df.append(df_dist)
    x_first =  df_tot[df_tot['Fold'] == args.pdb1]
    x_first =  x_first[['j','i']].to_numpy()
    x_second = df_tot[df_tot['Fold'] == args.pdb2]
    x_second = x_second[['i','j']].to_numpy()
    x_common = df_tot[df_tot['Fold'] == 'both']
    x_common = x_common[['j','i']].to_numpy()
    print('-------Fold-Switching crystalographic information-------')
    print('Number of unique contacts in PDB1:        {:}'.format(int(len(x_first))))
    print('Number of unique contacts in PDB2:        {:}'.format(int(len(x_second))))
    print('Number of unique contacts common to both: {:}\n'.format(int(len(x_common))))
    print('\n')
    
    df_sorted = db.Run(x_first,x_second,x_common,args.msa)
    df_sorted.to_csv('df_filtered.csv')
    
    param = db.EPS()
    print('DBSCAN parameters: eps = {:}  |  min = {:}'.format(param[0],param[1]))
    print('\n')
    plot = PLOT()
    plot.FS_plot(args.pdb1, args.pdb2, df_dist, df, args.msa, df_sorted)
    print('-------Results for full Pipeline-------')
    print('Total predictions: {:}'.format(len(df_sorted[df_sorted['group'] != -1].index)))
    print('Percent within +-1 of first fold:  {:}'.format(len(df_sorted[(df_sorted['group'] != -1) & (df_sorted['sort'] == 'pdb_1')])))
    print('Percent within +-1 of second fold: {:}'.format(len(df_sorted[(df_sorted['group'] != -1) & (df_sorted['sort'] == 'pdb_2')])))
    print('Percent within +-1 of common fold: {:}'.format(len(df_sorted[(df_sorted['group'] != -1) & (df_sorted['sort'] == 'common')])))
    print('Average radial distance of noise from xcrystal : {:}'.format(df_sorted['r'].mean()))
    print('\n')
    
    
    #___________________________________________________________________________________________________
    #   Generate 6-panel figure showing differences in GREMLIN algorithm and MSA-Transformer
    #   Also generates bar graph showing change in signal and noise and runs the hypergeometric test
    #___________________________________________________________________________________________________
   
    if args.Extra == 'y':
        extra = GMN_Extra()
        extra.MSA_Original_gmn(args.msa, xcontact, total_dist,args.Prediction_Offset,alignment,x_first,x_second,x_common)
        extra.Super_No_Filter(args.msa, xcontact, total_dist, args.Prediction_Offset, alignment, x_first, x_second, x_common,pdbs.N)
        extra.Subfamily_MSA_GMN(args.msa, xcontact, total_dist, args.Prediction_Offset, alignment,x_first,x_second,x_common)
        extra.MSA_Original_msatr(args.msa, xcontact, total_dist, args.Prediction_Offset, alignment, x_first, x_second, x_common)
        extra.Subfamily_MSA_MSATR(args.msa, xcontact, total_dist, args.Prediction_Offset, alignment, x_first, x_second, x_common)
        extra.GMN_MSATR(args.msa, xcontact, total_dist, args.Prediction_Offset, alignment, x_first, x_second, x_common, pdbs.N)
        extra.PLOT_EXTRA(args.msa, args.pdb1, args.pdb2,df_dist,df,df_sorted)
        gmn_file = Path('{:}/df_gmn_00.csv'.format(args.msa[:-4]))
        if gmn_file.is_file() :
            hp = HYPERGEO(args.pdb1, args.pdb2, args.msa, extra.DF_GMN_MSATR())
            hp.P_value(extra.DF_no_filter(),extra.DF_original(),extra.DF_orig_msatr(),pdbs.N, x_first, x_second, x_common)
