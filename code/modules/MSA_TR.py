#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code runs the MSATransformer algorithm on MSAs passed to it. This is adapted from https://github.com/rmrao/msa-transformer

@author: schaferjw
"""
import esm
import torch
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


def Make_DataFrame(tensor_in,output,query_length):
    mrf = tensor_in
    
    temp = {'i':[],'j':[],'zscore':[]}
    for i in range(mrf.shape[0]):
        for j in range(mrf.shape[1]):
            temp['i'].append(i)
            temp['j'].append(j)
            temp['zscore'].append(mrf[i][j].item())
      
    pd_mtx = pd.DataFrame.from_dict(temp)
    
    temp = pd_mtx
    N = int(query_length)
    zero = np.zeros((N,N))
    for index, row in temp.iterrows():
        zero[int(row['i']),int(row['j'])] = row['zscore']
        zero[int(row['j']),int(row['i'])] = row['zscore']
    
    top = pd_mtx.loc[pd_mtx['j'] - pd_mtx['i'] > 3].sort_values("zscore",ascending=False)

    return top.head(int(output*mrf.shape[0])),zero

def MSA_TR(msa,files,output,query_length,Beg=0,End=0):
    #run msa transformer on all alignments
    print('Start msa transformer...')
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
        top,full = Make_DataFrame(msa_contacts[n], output[n],query_length)
        top['i'] = top['i'] + Beg
        top['j'] = top['j'] + Beg
        np.savetxt('msatr/full_{:}.csv'.format(files[n][:-4]), full, delimiter=",")


            
