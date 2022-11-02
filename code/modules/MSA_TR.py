#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is adapted from https://github.com/rmrao/msa-transformer

@author: schaferjw
"""
import esm
import torch
import difflib
from Bio import SeqIO
import itertools
from typing import List, Tuple
import string
from scipy import stats
import numpy as np
import pandas as pd
torch.set_grad_enabled(False)

# This is an efficient way to delete lowercase characters and insertion characters from a string
deletekeys = dict.fromkeys(string.ascii_lowercase)
deletekeys["."] = None
deletekeys["*"] = None
translation = str.maketrans(deletekeys)

def read_sequence(filename: str) -> Tuple[str, str]:
    """ Reads the first (reference) sequences from a fasta or MSA file."""
    record = next(SeqIO.parse(filename, "fasta"))
    return record.description, str(record.seq)

def remove_insertions(sequence: str) -> str:
    """ Removes any insertions into the sequence. Needed to load aligned sequences in an MSA. """
    return sequence.translate(translation)

def read_msa(filename: str, nseq: int) -> List[Tuple[str, str]]:
    """ Reads the first nseq sequences from an MSA file, automatically removes insertions."""
    return [(record.description, remove_insertions(str(record.seq)))
            for record in itertools.islice(SeqIO.parse(filename, "fasta"), nseq)]

def normalize(x):
  x = stats.boxcox(x - np.amin(x) + 1.0)[0]
  x_mean = np.mean(x)
  x_std = np.std(x)
  return((x-x_mean)/x_std)


def Make_DataFrame(tensor_in,output):
    mrf = tensor_in
    
    temp = {'i':[],'j':[],'zscore':[]}
    for i in range(mrf.shape[0]):
        for j in range(mrf.shape[1]):
            temp['i'].append(i)
            temp['j'].append(j)
            temp['zscore'].append(mrf[i][j].item())
      
    pd_mtx = pd.DataFrame.from_dict(temp)
    
    top = pd_mtx.loc[pd_mtx['j'] - pd_mtx['i'] > 3].sort_values("zscore",ascending=False)

    return top.head(int(output*mrf.shape[0]))

def MSA_TR(msa,files,Beg=0,End=0):
    #run msa transformer on all alignments
    print('Start msa transformer...')
    output = np.array([])
    for i in files:
        num_temp = sum(1 for line in open('{:}/{:}'.format(msa,i)))
        output = np.append(output,num_temp)
    output = ((2.-1.5)/(np.max(output)-np.min(output)))*(output-np.min(output)) + 1.5
    if np.isnan(output[0]) == True:
        output[0] = 2.0
    msa_data = []
    for i in files:
        temp = read_msa("{:}/{:}".format(msa,i),64)
        if Beg != 0:
            temp_edit = []
            for i in range(len(temp)):
                short = (temp[i][0],temp[i][1][Beg:])
                temp_edit.append(short)
            temp = temp_edit
        if End != 0:
            temp_edit = []
            for i in range(len(temp)):
                short = (temp[i][0],temp[i][1][:-End])
                temp_edit.append(short)
            temp = temp_edit
        msa_data.append(temp)
    
    msa_transformer, msa_alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    msa_batch_converter = msa_alphabet.get_batch_converter()
    msa_batch_labels, msa_batch_strs, msa_batch_tokens = msa_batch_converter(msa_data)
    print(msa_batch_tokens.size(), msa_batch_tokens.dtype)  # Should be a 3D tensor with dtype torch.int64.

    msa_contacts = msa_transformer.predict_contacts(msa_batch_tokens).cpu()
    
    for n in range(len(files)):
        top = Make_DataFrame(msa_contacts[n], output[n])
        top['i'] = top['i'] + Beg
        top['j'] = top['j'] + Beg
        top.to_csv('{:}/df_msatr_{:}.csv'.format(msa,files[n][:-4]))

            
