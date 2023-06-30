#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:58:58 2022

@author: schaferjw

Create a dual-fold contact map for a fold-switching protein.
Two structures are needed to run this, the structures need to have indices that are similar (i.e. if one chain goes from 25-200 and the other goes from 1000-1075 this will fail)
This script will find an optimal alignment of the contact maps between the two structures.

"""

from collections import Counter
import mdtraj as md
import numpy as np
import pandas as pd
import re

class PDB_STRUCTURE():
    def __init__(self,structure):
        print('\n')
        print('Loading crystallographic information.')
        # This section loads in the two pdbs and saves information needed for identifying intra-chain contacts
        self.pdb = md.load_pdb('{:}'.format(structure))
        self.topo = self.pdb.topology
        table,bonds = self.topo.to_dataframe()
        self.indices = np.array(self.topo.select('protein and chainid 0')) #index of atom number
        self.pdb = self.pdb.atom_slice(self.indices)
        self.last_res = re.search('[0-9]+', str(self.topo.atom(self.indices[-1]))) #identify last residue 
        self.first_res = re.search('[0-9]+', str(self.topo.atom(self.indices[1]))) #identify first residue
        
        #_______________________________________________________________________________________________________________________________

        self.N = len(table['resSeq'].unique())

        self.structure = structure

    def N(self):return self.N

    def Get_contacts(self,b_factor,distance):
        pairs,bfactor = [],[]
        #read every contact that is within distance (set to 8A by default) 
        #this code finds all contacts by atom identifies residue numbers and removes redunant contacts from list
        #this code also removes any contacts within 3 residues of origin
        for atom in zip(self.indices,self.pdb.xyz[0,self.indices,2]):
            contacts = md.compute_neighbors(self.pdb,cutoff=distance,query_indices=[atom[0]]) #contacts within 5A of index
            origin = re.search('[0-9]+',str(self.topo.atom(atom[0])))
            b = atom[1]
            for i in zip(contacts[0],self.pdb.xyz[0,contacts,2][0]):
                name = re.search('[0-9]+', str(self.topo.atom(i[0]))) #match the id of the sequence
                if int(name.group()) < int(origin.group())-2 or int(name.group()) > int(origin.group())+2:
                    if int(name.group()) > int(origin.group()):
                        if not [int(origin.group()),int(name.group())] in pairs:
                            pairs.append([int(origin.group()),int(name.group())])
                            bfactor.append(min(b,i[1]))
        
        #_______________________________________________________________________________________________________________________________

        if b_factor == 'y':
            pairs = [[contact[0][0],contact[0][1],contact[1]] for contact in zip(pairs,bfactor)]
        return pairs
        #_______________________________________________________________________________________________________________________________

    def Multimer_contact(self,b_factor,distance):
        # This section loads in the two pdbs and saves information needed for identifying inter-chain contacts
        multi = md.load_pdb('{:}'.format(self.structure))
        multi_topo = multi.topology
        multi_indices = np.array(multi_topo.select('protein and chainid 0')) #index of atom number
        last_res = re.search('[0-9]+', str(multi_topo.atom(multi_indices[-1]))) #identify last atom index number
        pairs = []
        chain_pairs = []
        
        #_______________________________________________________________________________________________________________________________

        #read every contact that is within distance (set to 10A by default) 
        #this code finds all contacts by atom identifies residue numbers and removes redunant contacts from list
        #this code also removes any contacts within 3 residues of origin            
        for atom in zip(multi_indices):
            contacts = md.compute_neighbors(multi,cutoff=distance,query_indices=[atom[0]]) #contacts within 10A of index
            origin = re.search('[0-9]+',str(multi_topo.atom(atom[0])))
            for i in zip(contacts[0]):
                name = re.search('[0-9]+', str(multi_topo.atom(i[0]))) #match the id of the sequence
                if int(name.group()) < int(origin.group())-4 or int(name.group()) > int(origin.group())+4 and int(name.group()) <int(last_res.group())+1:
                    if int(name.group()) > int(origin.group()):
                        if not [int(origin.group()),int(name.group())] in pairs:
                            pairs.append([int(origin.group()),int(name.group())])
                if i > multi_indices[-1]: #identify any residue number larger than origin chains largest index
                    if int(name.group()) > int(origin.group()):
                        if not (int(origin.group()),int(name.group())) in chain_pairs and not (int(origin.group()),int(name.group())) in pairs:
                            if len(name.group()) >3:
                                chain_pairs.append((int(origin.group()),int(name.group()[1:])))
                            else:
                                chain_pairs.append((int(origin.group()),int(name.group())))

        #_______________________________________________________________________________________________________________________________
        #create a separate  list of contacts within 10A that are between chain A and another chain (this will be empty in pdb's that only have one chain)
        chain_pairs = list(set(chain_pairs)) #using set is an efficient way of dealing with duplicates, note it does not retain the order of the original list
        chain_pairs = [self.reverse(i) if self.triangle(i[1],i[0]) > 0  else i for i in chain_pairs]
        chain_pairs = [[i[0],i[1]] for i in chain_pairs if i[1]-3 > i[0]]
        chain_pairs = chain_pairs #+ pairs  #pairs are intrachain contacts     
   
        return chain_pairs
        #_______________________________________________________________________________________________________________________________

            
    def triangle(self,x,y):return y-x
    def reverse(self,x):return x[::-1]
    
       
class xOPT():
    def __init__(self,contacts_pdb1,contacts_pdb2,pdb1N,pdb2N,structure_1,structure_2):
   	#search for optimal alignment of two pdb structures contacts
        print('Searching for best alignment of PDBs...')
        opt = [0,0]
        contacts_pdb1 = [[contact[0],contact[1]] for contact in contacts_pdb1]
        contacts_pdb2 = [[contact[0],contact[1]] for contact in contacts_pdb2]
        
        if pdb1N < pdb2N:
            for increment in range(-int(pdb1N/2),int(pdb1N/2)):
                
                a,b = [],[]
                pairs = [[i[0]+increment,i[1]+increment] for i in contacts_pdb1]       
                for i in pairs:
                   if i in contacts_pdb2:a.append('common')    #contact exists in both folds
                   else:a.append('{:}'.format(structure_1)) #contact only exists in fold 1
    
                for i in contacts_pdb2:
                   if i in pairs:b.append('common')          #contact exists in both folds
                   else:b.append('{:}'.format(structure_2)) #contact only exists in fold 2
                concat = a + b       
                count = Counter(concat) #creates a dictionary of each designation with the associated count
                if opt[0] < count['common']:
                    opt = (count['common'], increment) #optimal alignment is stored here as (number_of_both, number_to_add_to_index)
        else:
            for increment in range(-int(pdb2N/2),int(pdb2N/2)):
                a,b = [],[]
                pairs = [[i[0]+increment,i[1]+increment] for i in contacts_pdb2]
                for i in pairs:
                   if i in contacts_pdb1:a.append('common')         #contact exists in both folds
                   else:a.append('{:}'.format(structure_2)) #contact only exists in fold 1
    
                for i in contacts_pdb1:
                   if i in pairs:b.append('common')     #contact exists in both folds
                   else:b.append('{:}'.format(structure_1)) #contact only exists in fold 2
                concat = a + b 
                count = Counter(concat) #creates a dictionary of each designation with the associated count
                if opt[0] < count['common']:
                    opt = (count['common'], increment) #optimal alignment is stored here as (number_of_both, number_to_add_to_index)
        
        print('Best alignment between pdbs found with adjustment of {:} yielding {:} redundant contacts.\n'.format(opt[1],opt[0]))

        self.opt = opt

    def OPT(self,manual,msa,contacts_pdb1,contacts_pdb2,pdb1N,pdb2N,structure_1,structure_2,b_factor):
       if b_factor == 'y':
            plddt_pdb1 = [contact[2] for contact in contacts_pdb1]
            plddt_pdb2 = [contact[2] for contact in contacts_pdb2]
       
       a,b = [],[]
       #add the perturbation to create the optimal alignment to the shorter pdb
       if manual == 'n':
           if pdb1N < pdb2N:
               contacts_pdb1 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in contacts_pdb1]
               contacts_pdb2 = [(i[0],i[1]) for i in contacts_pdb2]
           else:
               contacts_pdb2 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in contacts_pdb2]
               contacts_pdb1 = [(i[0],i[1]) for i in contacts_pdb1]

       else:
           if pdb1N < pdb2N:
               contacts_pdb1 = [(i[0]+int(manual),i[1]+int(manual)) for i in contacts_pdb1]
               contacts_pdb2 = [(i[0],i[1]) for i in contacts_pdb2]
           else:
               contacts_pdb2 = [(i[0]+int(manual),i[1]+int(manual)) for i in contacts_pdb2]
               contacts_pdb1 = [(i[0],i[1]) for i in contacts_pdb1]

       #_______________________________________________________________________________________________________________________________
       #recalculate the number of unique/nonunique contacts 
       for i in contacts_pdb1:
           if i in contacts_pdb2:a.append('both')  #this means this contact exitst in both folds
           else:a.append('{:}'.format(structure_1))  #this means this contact only exits in fold 1
       for i in contacts_pdb2:
           if i in contacts_pdb1:b.append('both')  #this means this contact exitst in both folds
           else:b.append('{:}'.format(structure_2))  #this means this contact only exits in fold 2
       #_______________________________________________________________________________________________________________________________
       if manual == 'n':
           if pdb1N < pdb2N:
               contacts_pdb2 = [(i[1],i[0]) for i in contacts_pdb2]
           else:
               contacts_pdb1 = [(i[1],i[0]) for i in contacts_pdb1]
       else:
           if pdb1N < pdb2N:
               contacts_pdb2 = [(i[1],i[0]) for i in contacts_pdb2]
           else:
               contacts_pdb1 = [(i[1],i[0]) for i in contacts_pdb1]
       concat = a + b
       
       # save xcrystal information
       self.switch = 'n'
       if contacts_pdb1[0][0] < contacts_pdb1[0][1]:
           contacts_pdb1 = [(contact[1],contact[0]) for contact in contacts_pdb1]
           contacts_pdb2 = [(contact[1],contact[0]) for contact in contacts_pdb2]
           self.switch = 'y'
       
       if b_factor == 'y':
           plddt = plddt_pdb1 + plddt_pdb2
       else:
           plddt = ['none' for i in range(len(concat))]
       xcontact = contacts_pdb1 + contacts_pdb2
       df = pd.DataFrame(xcontact,columns=list('ij'))
       df['Fold'] = concat
       df.to_csv('{:}/df_xcontact.csv'.format(msa))                
           
       return xcontact, self.opt[1], df, plddt
   
    def OPT_Multi(self,dist_contacts_pdb1,dist_contacts_pdb2,manual,msa,pdb1N,pdb2N,structure_1,structure_2,b_factor):
        if b_factor == 'y':
             plddt_pdb1 = [contact[2] for contact in dist_contacts_pdb1]
             plddt_pdb2 = [contact[2] for contact in dist_contacts_pdb2]
        a,b = [],[]
       #add the perturbation to create the optimal alignment to the shorter pdb
        if manual == 'n':
           if pdb1N < pdb2N:
               dist_contacts_pdb1 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in dist_contacts_pdb1]
               dist_contacts_pdb2 = [(i[0],i[1]) for i in dist_contacts_pdb2]
           else:
               dist_contacts_pdb2 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in dist_contacts_pdb2]
               dist_contacts_pdb1 = [(i[0],i[1]) for i in dist_contacts_pdb1]
        else:
           if pdb1N < pdb2N:
               dist_contacts_pdb1 = [(i[0]+int(manual),i[1]+int(manual)) for i in dist_contacts_pdb1]
               dist_contacts_pdb2 = [(i[0],i[1]) for i in dist_contacts_pdb2]
           else:
               dist_contacts_pdb2 = [(i[0]+int(manual),i[1]+int(manual)) for i in dist_contacts_pdb2]
               dist_contacts_pdb1 = [(i[0],i[1]) for i in dist_contacts_pdb1]

       #_______________________________________________________________________________________________________________________________
       #recalculate the number of unique/nonunique contacts 
        for i in dist_contacts_pdb1:
           if i in dist_contacts_pdb2:a.append('both')  #this means this contact exitst in both folds
           else:a.append('{:}'.format(structure_1))  #this means this contact only exits in fold 1
        for i in dist_contacts_pdb2:
           if i in dist_contacts_pdb1:b.append('both')  #this means this contact exitst in both folds
           else:b.append('{:}'.format(structure_2))  #this means this contact only exits in fold 2
       #_______________________________________________________________________________________________________________________________
        if manual == 'n':
           if pdb1N < pdb2N:
               dist_contacts_pdb2 = [(i[1],i[0]) for i in dist_contacts_pdb2]
           else:
               dist_contacts_pdb1 = [(i[1],i[0]) for i in dist_contacts_pdb1]
        else:
           if pdb1N < pdb2N:
               dist_contacts_pdb2 = [(i[1],i[0]) for i in dist_contacts_pdb2]
           else:
               dist_contacts_pdb1 = [(i[1],i[0]) for i in dist_contacts_pdb1]
        concat = a + b
        
        # save xcrystal information
        # if dist_contacts_pdb1[0][0] < dist_contacts_pdb1[0][1]:
        if self.switch == 'y':
            dist_contacts_pdb1 = [(contact[1],contact[0]) for contact in dist_contacts_pdb1]
            dist_contacts_pdb2 = [(contact[1],contact[0]) for contact in dist_contacts_pdb2]
       
        if b_factor == 'y':
            plddt = plddt_pdb1 + plddt_pdb2
        else:
            plddt = ['none' for i in range(len(concat))]
        
        # save xcrystal information 
        xcontact = dist_contacts_pdb1 + dist_contacts_pdb2
        df = pd.DataFrame(xcontact,columns=list('ij'))
        df['Fold'] = concat
        df.to_csv('{:}/df_dist_xcontact.csv'.format(msa))
        return xcontact, df
    