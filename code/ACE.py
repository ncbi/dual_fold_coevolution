#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:32:20 2022

@author: schaferjw

This code is intended for evaluating the overlap between crystallographic contacts of known fold-switching proteins
and coevolutionary predictions made using GREMLIN and MSA Transformer.

To run this code .pdb files are required for each conformation and an MSA file is needed in stockholm format.
NOTE: the stockholm files were generated with HMMER 3.3.2 and renamed as 'pdbID.msa' (this program expects the .msa extension
"""

import argparse, os
import numpy as np
from Bio import AlignIO
from pathlib import Path
import pandas as pd

from modules.MSA_PREP import Edit_MSA,TARFILE
from modules.MSA_TR   import MSA_TR
from modules.SUPER import Data
from modules.COEVOLUTION import Start
from modules.PDB2 import PDB_STRUCTURE, xOPT
from modules.DB_SORT import DB_SCAN
from modules.PLOT import PLOT
from modules.EXTRA import GMN_Extra
from modules.HYPERGEO import HYPERGEO
from modules.SUPER_DECOMP import Decomp


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--msa", type=str,default='None', help='MSA file in stockholm format')
    parser.add_argument("--Sub_Family_Parser",type=str,default='QID',help='Choose between parsing sub-family alignments by QID or E-value (E-value is from HMMER)')
    parser.add_argument("--fasta", type=str,default='None', help='MSA file in fasta format')
    parser.add_argument("--b_factor",type=str,default='n',help="(y/n) ESMfold and Alphafold have pLDDT scores in the b_factor column. This flag colors the dualfold contact map based on pLDDT, poor predictions are red and ok predictions are yellow.")
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
    
    #___________________________________________________________________________________________________  
    #   Prepare Subfamily MSAs from alignment generated from HMMER
    #___________________________________________________________________________________________________
    
    if args.fasta == 'None':pass
    elif args.msa == 'None':args.msa = args.fasta
    else:
        print('An MSA must be provided to run coevolutionary analysis...')
        exit()

    if os.path.isdir(args.msa[:-4]) == False:        
        os.system('mkdir {:}'.format(args.msa[:-4]))
    if args.fasta == 'None':align = AlignIO.read("{:}".format(args.msa), "stockholm")
    else:align = AlignIO.read("{:}".format(args.msa), "fasta")
    
    seqs = [str(record.seq[:]) for record in align]
    query_length = len(seqs[0].replace('-',''))
    names = [record.id[:35].ljust(36) for record in align]

    if os.path.isdir('tmp') == True: 
        TARFILE.Open_File(args.msa[:-4],'msatr.tar.gz')
        TARFILE.Open_File(args.msa[:-4],'gmn.tar.gz')
        TARFILE.Open_File(args.msa[:-4],'csv.tar.gz')
    else:
        edit_msa = Edit_MSA(args.msa,seqs,names,args.Sub_Family_Parser)


    #___________________________________________________________________________________________________
    #   Run msa transformer
    #___________________________________________________________________________________________________
    
    files = [i for i in os.listdir('{:}/'.format(args.msa[:-4])) if os.path.isfile(os.path.join('{:}/'.format(args.msa[:-4]), i))]
    files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])
    remove = []
    files = [file for file in files if file[-4:] == '.a3m']
    for i in files:
        num_temp = sum(1 for line in open('{:}/{:}'.format(args.msa[:-4],i)))
        if num_temp/2 < 100 and i[0:2] != '00': 
            os.system('rm -f {:}/{:}'.format(args.msa[:-4],i))
            remove.append(i)
    files_tr = [file for file in files if file not in remove]
    output = [sum(1 for line in open('{:}/{:}'.format(args.msa[:-4],file))) for file in files_tr]
    if output: 
        output = ((2.-1.5)/(np.max(output)-np.min(output)))*(output-np.min(output)) + 1.5
        if np.isnan(output[0]) == True:output[0] = 2.0    
    if os.path.isdir('msatr') == False: 
        os.system('mkdir {:}'.format('msatr'))
        MSA_TR(args.msa[:-4],files_tr,output,query_length,args.Beg,args.End)
    else:
        files_complete = [i for i in os.listdir('msatr') if os.path.isfile(os.path.join('msatr', i))]
        files_complete = sorted(files_complete, key = lambda x: x.rsplit('.', 1)[0])
        files_complete = [file for file in files_tr if any(complete[5:11] in file for complete in files_complete)]
        files_msatr_idx = [idx for idx in range(len(files_tr)) if files_tr[idx] not in files_complete]
        files_tr = [file for file in files_tr if file not in files_complete]
        output = [output[idx] for idx in files_msatr_idx]
        if files_tr:MSA_TR(args.msa[:-4],files_tr,output,query_length,args.Beg,args.End) 
    
    
    #___________________________________________________________________________________________________
    #   Run Coevolutionary analysis
    #___________________________________________________________________________________________________
    
    files = [i for i in os.listdir('{:}/'.format(args.msa[:-4])) if os.path.isfile(os.path.join('{:}/'.format(args.msa[:-4]), i))]
    files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])
    files_gmn = [file for file in files if file[-4:] == '.out']
    if files_gmn:
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
                remove.append(i)
        files_gmn = [file for file in files_gmn if file not in remove]
        remove = []
        if not files_gmn:
            print('...Not enough data for GREMLIN algorithm...')
        else:
            output = [sum(1 for line in open('{:}/{:}'.format(args.msa[:-4],file))) for file in files_gmn]
            if os.path.isdir('gmn') == False: 
                os.system('mkdir {:}'.format('gmn'))
                output = ((args.Depth_max-args.Depth_min)/(np.max(output)-np.min(output)))*(output-np.min(output)) + args.Depth_min
                if np.isnan(output[0]) == True:
                    output[0] = 2.0
                for file in zip(files_gmn,np.reshape(output,(len(files_gmn),1)).tolist()):
                    print('Coevolutionary Analysis running on {:}'.format(file[0]))
                    full = Start(args.msa[:-4],file[0],file[1][0],query_length)
                    np.savetxt('gmn/full_{:}.csv'.format(file[0][:-4]), full, delimiter=",")
            else:
                output = ((args.Depth_max-args.Depth_min)/(np.max(output)-np.min(output)))*(output-np.min(output)) + args.Depth_min
                if np.isnan(output[0]) == True:
                    output[0] = 2.0
                files_complete = [i for i in os.listdir('gmn') if os.path.isfile(os.path.join('gmn', i))]
                files_complete = sorted(files_complete, key = lambda x: x.rsplit('.', 1)[0])
                files_complete = [file for file in files_gmn if any(complete[5:11] in file for complete in files_complete)]
                files_gmn_idx = [idx for idx in range(len(files_gmn)) if files_gmn[idx] not in files_complete]
                files_gmn = [file for file in files_gmn if file not in files_complete]
                output = [output[idx] for idx in files_gmn_idx]
                if files_gmn:
                    for file in zip(files_gmn,np.reshape(output,(len(files_gmn),1)).tolist()):
                        print('Coevolutionary Analysis running on {:}'.format(file[0]))
                        full = Start(args.msa[:-4],file[0],file[1][0],query_length)
                        np.savetxt('gmn/full_{:}.csv'.format(file[0][:-4]), full, delimiter=",")


    elif not files_gmn:
            print('...Not enough data for GREMLIN algorithm...')
     
       
   
    #___________________________________________________________________________________________________
    #   Superimpose predictions into one contact map
    #___________________________________________________________________________________________________
    
    files = [i for i in os.listdir('{:}/'.format(args.msa[:-4])) if os.path.isfile(os.path.join('{:}/'.format(args.msa[:-4]), i))]
    files = sorted(files, key = lambda x: x.rsplit('.', 1)[0])
    msas_msatr = [file for file in files if file[-4:] == '.a3m']
    msas_gmn =   [file for file in files if file[-4:] == '.out']
    files_csv =   [file for file in files if file[-4:] == '.csv' and file[:2] != 'df']

    files_msatr2 = [i for i in os.listdir('msatr') if os.path.isfile(os.path.join('msatr', i))]
    files_msatr2 = sorted(files_msatr2, key = lambda x: x.rsplit('.', 1)[0])
    files_msatr2 = [file for file in files_msatr2 if file[-4:] == '.csv']
    files_msatr2 = ['msatr/'+file for file in files_msatr2]

    files_gmn2 = [i for i in os.listdir('gmn') if os.path.isfile(os.path.join('gmn', i))]
    files_gmn2 = sorted(files_gmn2, key = lambda x: x.rsplit('.', 1)[0])
    files_gmn2 = [file for file in files_gmn2 if file[-4:] == '.csv']
    files_gmn2 = ['gmn/'+file for file in files_gmn2]   
    
    files_all = files_gmn2 + files_msatr2
    predictions = Data(files_all)
    df_super = predictions.Most_Probable()
    
    #___________________________________________________________________________________________________
    #   Load crystallagraphic information
    #___________________________________________________________________________________________________
    #MDtraj uses deprecated np.int instead of the new np.int32 or np.int64
    #comment these lines out to see full deprecation warning
    import warnings
    warnings.filterwarnings('ignore')
    
    pdb1 = PDB_STRUCTURE(args.pdb1)
    contacts_pdb1 = pdb1.Get_contacts(args.b_factor,0.8)       #define list of intrachain contacts for comparison
    dist_contacts_pdb1 = pdb1.Multimer_contact(args.b_factor,1.0)   #create list of inerchain contacts (based on above)
    pdb2 = PDB_STRUCTURE(args.pdb2)
    contacts_pdb2 = pdb2.Get_contacts(args.b_factor,0.8)       #define list of intrachain contacts for comparison
    dist_contacts_pdb2 = pdb2.Multimer_contact(args.b_factor,1.0)   #create list of inerchain contacts (based on above)
    
    xopt = xOPT(contacts_pdb1,contacts_pdb2,pdb1.N,pdb2.N,args.pdb1,args.pdb2)
    xcontact, opt, df_intra, plddt = xopt.OPT(args.Manual_PDB,args.msa[:-4],contacts_pdb1,contacts_pdb2,pdb1.N,pdb2.N,args.pdb1,args.pdb2,args.b_factor)
    dist_xcontact, df_dist = xopt.OPT_Multi(dist_contacts_pdb1,dist_contacts_pdb2,args.Manual_PDB,args.msa[:-4],pdb1.N,pdb2.N,args.pdb1,args.pdb2,args.b_factor)
    #___________________________________________________________________________________________________
    #   Find Clusters within predictions
    #   align superposition of predictions to crystallographic information and cluster
    #   Print information about xcrystal alignment, predictions, and save png file.
    #___________________________________________________________________________________________________ 
   
    df_tot = df_intra.append(df_dist)
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
   
    df_super = predictions.Most_Probable(N_r=7.5)
    df_super.to_csv(f'{args.msa[:-4]}/df_super.csv')
    if args.Prediction_Offset != 'n':
        df_super['i'] = df_super['i'] + int(args.Prediction_Offset)
        df_super['j'] = df_super['j'] + int(args.Prediction_Offset)
    
    db = DB_SCAN(df_super,xcontact,args.Manual_Pred)
    alignment = db.Return_Alignment()
    df_sorted = db.Run(x_first,x_second,x_common,args.msa,query_length)
    df_sorted.to_csv('df_filtered.csv')
    param = db.EPS()
    print('DBSCAN parameters: eps = {:}  |  min = {:}'.format(param[0],param[1]))
    print('\n')
    plot = PLOT()
    plot.FS_plot(args.pdb1, args.pdb2, df_dist, df_intra, args.msa, df_sorted,plddt=plddt)
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
        extra.MSA_Original_gmn(args.msa, xcontact, dist_xcontact,args.Prediction_Offset,alignment,x_first,x_second,x_common,query_length)
        extra.Super_No_Filter(args.msa, xcontact, dist_xcontact, args.Prediction_Offset, alignment, x_first, x_second, x_common,query_length)
        extra.Subfamily_MSA_GMN(args.msa, xcontact, dist_xcontact, args.Prediction_Offset, alignment,x_first,x_second,x_common,query_length)
        extra.MSA_Original_msatr(args.msa, xcontact, dist_xcontact, args.Prediction_Offset, alignment, x_first, x_second, x_common,query_length)
        extra.Subfamily_MSA_MSATR(args.msa, xcontact, dist_xcontact, args.Prediction_Offset, alignment, x_first, x_second, x_common,query_length)
        extra.GMN_MSATR(args.msa, xcontact, dist_xcontact, args.Prediction_Offset, alignment, x_first, x_second, x_common)
        extra.PLOT_EXTRA(args.msa, args.pdb1, args.pdb2,df_dist,df_intra,df_sorted,plddt=plddt)
        gmn_file = Path('{:}/df_gmn_00.csv'.format(args.msa[:-4]))
        hp = HYPERGEO(args.pdb1, args.pdb2, args.msa, extra.DF_GMN_MSATR())
        hp.P_value(extra.DF_no_filter(),extra.DF_original(),extra.DF_orig_msatr(),max(pdb1.N,pdb2.N), x_first, x_second, x_common)

    #___________________________________________________________________________________________________
    #   subfamilies by iteratively removing larger alignments
    #___________________________________________________________________________________________________
    

        gmn_decomp = Decomp(files_gmn2,'GMN')
        gmn_decomp.Pred(args.msa,xcontact,msas_gmn,df_intra,df_dist,args.pdb1,args.pdb2,args.Prediction_Offset,alignment,'gmn',query_length)
    
        msatr_decomp = Decomp(files_msatr2,'MSATR')
        msatr_decomp.Pred(args.msa,xcontact,msas_msatr,df_intra,df_dist,args.pdb1,args.pdb2,args.Prediction_Offset,alignment,'msatr',query_length)

        files_full = []
        #need a list of alternating gmn and msatr...
        for idx in range(max(len(files_gmn2),len(files_msatr2))):
            try:files_full.append(files_gmn2[idx])
            except:print('')
            try:files_full.append(files_msatr2[idx])
            except:print('')
        full_decomp = Decomp(files_full,'FULL')
        full_decomp.Pred(args.msa,xcontact,msas_msatr,df_intra,df_dist,args.pdb1,args.pdb2,args.Prediction_Offset,alignment,'full',query_length)
 
    #___________________________________________________________________________________________________
    #   Compress intermediate files
    #___________________________________________________________________________________________________
    if os.path.isdir('tmp') == False: 
        os.system('mkdir tmp')
   
    TARFILE.Compress(args.msa[:-4],msas_msatr,'msatr')
    TARFILE.Compress(args.msa[:-4],msas_gmn,'gmn')
    TARFILE.Compress(args.msa[:-4],files_csv,'csv')