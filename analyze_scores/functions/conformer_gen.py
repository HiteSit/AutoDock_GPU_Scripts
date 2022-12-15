# Author: Dr. Riccardo Fusco
# Copyright Riccardo Fusco 2022

#!/usr/bin/env python
# coding: utf-8

import glob, os, shutil
import re
import pdb

import pandas as pd
import numpy as np

from openbabel import pybel, openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import chain

def gen_conf(smile, num_conf, num_iter):
    '''
    
    Generate and minimize some conformers, return a tuple


    '''
    mol = Chem.MolFromSmiles(smile)
    mol_h_MMFF = Chem.AddHs(mol)

    # Number of conformers to be generated
    num_of_conformer = num_conf
    max_iter = num_iter

    # Default value for min energy conformer
    min_energy_MMFF = 10000
    min_energy_index_MMFF = 0

    # Generate conformers (stored inside of the mol object)
    cids = AllChem.EmbedMultipleConfs(mol_h_MMFF, numConfs = num_of_conformer, params = AllChem.ETKDG())

    # Indexing conformes
    ids = list(cids)

    # Result minimized conformers (tuple)
    results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(mol_h_MMFF,max_iter,)

    # Find global minimium
    ndx = []
    energy = []
    
    for index, result in enumerate(results_MMFF):
        # if result[0] == 0:
        ndx.append(index)
        energy.append(result[1])
        
        # else:
        #     ndx.append(index)
        #     energy.append(np.nan)
    
    energy_round = [round(x,2) for x in energy]
    tup = list(zip(ndx,energy_round))
    
    min_conf = min(tup, key=lambda t: t[1])
    return min_conf

def convert_to_mol2(pdbqt_string):
    '''
    
    Convert string from pdbqt to mol2, return a mol2 string
    

    '''
    
    converted = []
    
    for run in pdbqt_string:
        mol = pybel.readstring("pdbqt",run)
        mol_string = mol.write("mol2")
        
        mol = pybel.readstring("mol2",mol_string)
        mol.addh()
        converted.append(mol.write("mol2"))
        
    return converted

def calc_conf_energy(mol2_string):
    '''
    
    Calculate conformational energy of a ligand pose, return a integer
    

    '''

    mol = Chem.MolFromMol2Block(mol2_string,removeHs=False,sanitize=False)
    mol_h = Chem.AddHs(mol)
    # mol_h_min = AllChem.MMFFOptimizeMolecule(mol_h)
    ff = AllChem.MMFFGetMoleculeForceField(mol_h, AllChem.MMFFGetMoleculeProperties(mol_h))
    ff.Initialize()
    ff.CalcEnergy()
    
    return ff.CalcEnergy()
