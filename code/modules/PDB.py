#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:58:58 2022

@author: schaferjw
"""

from collections import Counter
import mdtraj as md
import numpy as np
import pandas as pd
import re

class PDB_STRUCTURE():
    def __init__(self,structure_1,structure_2):
        print('\n')
        print('Loading crystallographic information.')
        # This section loads in the two pdbs and saves information needed for identifying intra-chain contacts
        self.pdb = md.load_pdb('{:}'.format(structure_1))
        self.topo = self.pdb.topology
        table,bonds = self.topo.to_dataframe()
        len_pdb1 = len(table['resSeq'].unique())
        self.indices = np.array(self.topo.select('protein and chainid 0')) #index of atom number
        self.pdb = self.pdb.atom_slice(self.indices)
        self.last_res = re.search('[0-9]+', str(self.topo.atom(self.indices[-1]))) #identify last residue 
        self.first_res = re.search('[0-9]+', str(self.topo.atom(self.indices[1]))) #identify first residue
        
        self.pdb2 = md.load_pdb('{:}'.format(structure_2)) 
        self.topo2 = self.pdb2.topology
        table,bonds = self.topo2.to_dataframe()
        len_pdb2 = len(table['resSeq'].unique())
        self.indices2 = np.array(self.topo2.select('protein and chainid 0')) #index of atom number
        self.pdb2 = self.pdb2.atom_slice(self.indices2)
        self.last_res2 = re.search('[0-9]+', str(self.topo2.atom(self.indices2[-1]))) #identify last residue
        self.first_res2 = re.search('[0-9]+', str(self.topo2.atom(self.indices2[1]))) #identify first residue
        #_______________________________________________________________________________________________________________________________

        self.N = max([len_pdb1,len_pdb2])

        self.structure_1 = structure_1
        self.structure_2 = structure_2

    def N(self):return self.N

    def Get_contacts(self,distance):
        pairs = []
        pairs2 = []
        #read every contact that is within distance (set to 8A by default) 
        #this code finds all contacts by atom identifies residue numbers and removes redunant contacts from list
        #this code also removes any contacts within 3 residues of origin
        for atom in zip(self.indices):
            contacts = md.compute_neighbors(self.pdb,cutoff=distance,query_indices=[atom[0]]) #contacts within 5A of index
            origin = re.search('[0-9]+',str(self.topo.atom(atom[0])))
            for i in zip(contacts[0]):
                name = re.search('[0-9]+', str(self.topo.atom(i[0]))) #match the id of the sequence
                if int(name.group()) < int(origin.group())-2 or int(name.group()) > int(origin.group())+2:
                    if not (int(origin.group()),int(name.group())) in pairs:
                        pairs.append((int(origin.group()),int(name.group())))

        for atom in zip(self.indices2):
            contacts2 = md.compute_neighbors(self.pdb2,cutoff=distance,query_indices=[atom[0]]) #contacts within 5A of index
            origin2 = re.search('[0-9]+',str(self.topo2.atom(atom[0])))
            for i in zip(contacts2[0]):
                name2 = re.search('[0-9]+', str(self.topo2.atom(i[0]))) #match the id of the sequence
                if int(name2.group()) < int(origin2.group())-2 or int(name2.group()) > int(origin2.group())+2:
                    if not (int(origin2.group()),int(name2.group())) in pairs2:
                        pairs2.append((int(origin2.group()),int(name2.group())))  
        
        #_______________________________________________________________________________________________________________________________
        # Seperate folds into upper and lower triangle
        pairs2 = [i for i in pairs2 if self.triangle(i[1],i[0]) < 0] #lower
        pairs =  [i for i in pairs  if self.triangle(i[1],i[0]) > 0] #upper

        #_______________________________________________________________________________________________________________________________

        self.pairs = pairs
        self.pairs2 = pairs2
            
    def triangle(self,x,y):return y-x
    def reverse(self,x):return x[::-1]
    
    def Combo_opt(self):
   	#search for optimal alignment of two pdb structures contacts
        print('Searching for best alignment of PDBs...')
        opt = [0,0]
        if self.pdb.n_residues < self.pdb2.n_residues:
            for j in range(-int(self.pdb.n_residues/3),int(self.pdb.n_residues/3)): 
                a,b = [],[]
                pairs = [(i[0]+j,i[1]+j) for i in self.pairs]                
                
                for i in pairs:
                   if i[::-1] in self.pairs2:a.append('both')    #contact exists in both folds
                   else:a.append('{:}'.format(self.structure_1)) #contact only exists in fold 1
    
                for i in self.pairs2:
                   if i[::-1] in pairs:b.append('both')          #contact exists in both folds
                   else:b.append('{:}'.format(self.structure_2)) #contact only exists in fold 2
                concat = a + b       
                count = Counter(concat) #creates a dictionary of each designation with the associated count
                if opt[0] < count['both']:
                    opt = (count['both'], j) #optimal alignment is stored here as (number_of_both, number_to_add_to_index)
        else:
            for j in range(-int(self.pdb2.n_residues/3),int(self.pdb2.n_residues/3)):
                a,b = [],[]
                pairs2 = [(i[0]+j,i[1]+j) for i in self.pairs2]

                for i in self.pairs:
                   if i[::-1] in pairs2:a.append('both')         #contact exists in both folds
                   else:a.append('{:}'.format(self.structure_1)) #contact only exists in fold 1
    
                for i in pairs2:
                   if i[::-1] in self.pairs:b.append('both')     #contact exists in both folds
                   else:b.append('{:}'.format(self.structure_2)) #contact only exists in fold 2
                concat = a + b       
                count = Counter(concat) #creates a dictionary of each designation with the associated count

                if opt[0] < count['both']:
                    opt = (count['both'], j) #optimal alignment is stored here as (number_of_both, number_to_add_to_index)
        print('Best alignment between pdbs found with adjustment of {:} yielding {:} redundant contacts.\n'.format(opt[1],opt[0]))

        self.opt = opt

    def Combo(self,manual,msa):
       a,b = [],[]
       #add the perturbation to create the optimal alignment to the shorter pdb
       if manual == 'n':
           if self.pdb.n_residues < self.pdb2.n_residues:
               self.pairs = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in self.pairs]
               self.dist_pdb_1 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in self.dist_pdb_1]
    
           else:
               self.pairs2 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in self.pairs2]
               self.dist_pdb_2 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in self.dist_pdb_2]
       else:
           if self.pdb.n_residues < self.pdb2.n_residues:
               self.pairs = [(i[0]+int(manual),i[1]+int(manual)) for i in self.pairs]
               self.dist_pdb_1 = [(i[0]+int(manual),i[1]+int(manual)) for i in self.dist_pdb_1]
    
           else:
               self.pairs2 = [(i[0]+int(manual),i[1]+int(manual)) for i in self.pairs2]
               self.dist_pdb_2 = [(i[0]+int(manual),i[1]+int(manual)) for i in self.dist_pdb_2]

       #_______________________________________________________________________________________________________________________________
       #recalculate the number of unique/nonunique contacts 
       for i in self.pairs:
           if i[::-1] in self.pairs2:
               a.append('both')  #this means this contact exitst in both folds
           else:
               a.append('{:}'.format(self.structure_1))  #this means this contact only exits in fold 1
       
       for i in self.pairs2:
           if i[::-1] in self.pairs:
               b.append('both')  #this means this contact exitst in both folds
           else:
               b.append('{:}'.format(self.structure_2))  #this means this contact only exits in fold 2
       #_______________________________________________________________________________________________________________________________
        
       concat = a + b
       c,d = [],[] 
       for i in self.dist_pdb_1:
           if i[::-1] in self.dist_pdb_2:
               c.append('both')
           else:
               c.append('{:}'.format(self.structure_1))
       for i in self.dist_pdb_2:
           if i[::-1] in self.dist_pdb_1:
               d.append('both')
           else:
               d.append('{:}'.format(self.structure_2))
       dist_concat = c + d
       # save xcrystal information 
       xcontact = self.pairs + self.pairs2
       df = pd.DataFrame(xcontact,columns=list('ij'))
       df['Fold'] = concat
       df.to_csv('{:}/df_xcontact.csv'.format(msa))
       total_dist = self.dist_pdb_1 + self.dist_pdb_2
       df_dist = pd.DataFrame(total_dist,columns=list('ij'))
       df_dist['Fold'] = dist_concat
       df_dist.to_csv('{:}/df_xcontact_dist.csv'.format(msa))
       return concat, dist_concat, self.opt[1], self.pairs, self.pairs2,self.dist_pdb_1,self.dist_pdb_2
       
    def Multimer_contact(self,distance):
        # This section loads in the two pdbs and saves information needed for identifying inter-chain contacts
        multi = md.load_pdb('{:}'.format(self.structure_1))
        multi_topo = multi.topology
        multi_indices = np.array(multi_topo.select('protein and chainid 0')) #index of atom number
        last_res = re.search('[0-9]+', str(multi_topo.atom(multi_indices[-1]))) #identify last atom index number
        pairs = []
        chain_pairs = []
        
        
        multi2 = md.load_pdb('{:}'.format(self.structure_2))
        multi2_topo = multi2.topology
        multi_indices2 = np.array(multi2_topo.select('protein and chainid 0')) #index of atom number
        last_res2 = re.search('[0-9]+', str(multi2_topo.atom(multi_indices2[-1]))) #identify last atom index number
        pairs2 = []
        chain_pairs2 = []
        #_______________________________________________________________________________________________________________________________

        #read every contact that is within distance (set to 10A by default) 
        #this code finds all contacts by atom identifies residue numbers and removes redunant contacts from list
        #this code also removes any contacts within 3 residues of origin            
        for atom in zip(multi_indices2):
            contacts2 = md.compute_neighbors(multi2,cutoff=distance,query_indices=[atom[0]]) #contacts within 10A of index
            origin2 = re.search('[0-9]+',str(multi2_topo.atom(atom[0])))
            for i in zip(contacts2[0]):
                name2 = re.search('[0-9]+', str(multi2_topo.atom(i[0]))) #match the id of the sequence
                if int(name2.group()) < int(origin2.group())-4 or int(name2.group()) > int(origin2.group())+4 and int(name2.group()) <int(last_res2.group())+1:
                    if not (int(origin2.group()),int(name2.group())) in pairs2:
                        pairs2.append((int(origin2.group()),int(name2.group())))
                if i > multi_indices2[-1]: #identify any residue number larger than origin chains largest index 
                        if not (int(origin2.group()),int(name2.group())) in chain_pairs2 and not (int(origin2.group()),int(name2.group())) in pairs2:
                            if len(name2.group()) >3:
                                chain_pairs2.append((int(origin2.group()),int(name2.group()[1:])))
                            else:
                                chain_pairs2.append((int(origin2.group()),int(name2.group())))

        for atom in zip(multi_indices):
            contacts = md.compute_neighbors(multi,cutoff=distance,query_indices=[atom[0]]) #contacts within 10A of index
            origin = re.search('[0-9]+',str(multi_topo.atom(atom[0])))
            for i in zip(contacts[0]):
                name = re.search('[0-9]+', str(multi_topo.atom(i[0]))) #match the id of the sequence
                if int(name.group()) < int(origin.group())-4 or int(name.group()) > int(origin.group())+4 and int(name.group()) <int(last_res.group())+1:
                    if not (int(origin.group()),int(name.group())) in pairs:
                        pairs.append((int(origin.group()),int(name.group())))
                if i > multi_indices[-1]: #identify any residue number larger than origin chains largest index 
                        if not (int(origin.group()),int(name.group())) in chain_pairs and not (int(origin.group()),int(name.group())) in pairs:
                            if len(name.group()) >3:
                                chain_pairs.append((int(origin.group()),int(name.group()[1:])))
                            else:
                                chain_pairs.append((int(origin.group()),int(name.group())))

        #_______________________________________________________________________________________________________________________________
        
        # make pairs an upper triangle and remove any contacts that are in intrachain list
        pairs =  [i for i in pairs  if self.triangle(i[1],i[0]) > 0 and i not in self.pairs] #upper

        #_______________________________________________________________________________________________________________________________
        #create a separate  list of contacts within 10A that are between chain A and another chain (this will be empty in pdb's that only have one chain)
        chain_pairs = list(set(chain_pairs)) #using set is an efficient way of dealing with duplicates, note it does not retain the order of the original list
        chain_pairs = [i if self.triangle(i[1],i[0]) > 0  else self.reverse(i) for i in chain_pairs]
        chain_pairs = [i for i in chain_pairs if i[0]-3 > i[1]]

        #_______________________________________________________________________________________________________________________________
        # make pairs a lower triangle and remove any contacts that are listed as being within 5A
        pairs2 =  [i for i in pairs2  if self.triangle(i[1],i[0]) < 0 and i not in self.pairs2] #upper

        #_______________________________________________________________________________________________________________________________

        #create a separate  list of contacts within 10A that are between chain A and another chain (this will be empty in pdb's that only have one chain)
        chain_pairs2 = list(set(chain_pairs2)) #using set is an efficient way of dealing with duplicates, note it does not retain the order of the original list
        chain_pairs2 = [i if self.triangle(i[1],i[0]) < 0  else self.reverse(i) for i in chain_pairs2 ]
        chain_pairs2 = [i for i in chain_pairs2 if i[0]+3 < i[1]]
        
        #_______________________________________________________________________________________________________________________________

        chain_pairs = chain_pairs + pairs
        chain_pairs2 = chain_pairs2 + pairs2
        
        
        if self.structure_1 == self.structure_2:
            chain_pairs2 = [self.reverse(i) for i in chain_pairs2]
            chain_pairs = list(set(chain_pairs+chain_pairs2))
            chain_pairs2 = [self.reverse(i) for i in chain_pairs]
                
        self.dist_pdb_1, self.dist_pdb_2 = chain_pairs,chain_pairs2
        #_______________________________________________________________________________________________________________________________
