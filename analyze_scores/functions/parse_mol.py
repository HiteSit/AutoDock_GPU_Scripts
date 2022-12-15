import glob, os, shutil
import re
import pdb

import pandas as pd
import numpy as np

from openbabel import pybel, openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import chain

def grab_mol_clean(file):
    '''
    
    Grab the numbers of the best clusters, find the right run, clean the "DOCKED:".
    So It write all the best runs in multiple pdbqt files.
    
    '''
        
    basename = file.split(".")[0]
    
    pattern_1 = r'[^\d][\s\d][\s\d]\s\W\s\s\s\s[\s\W][\s\d\W]\d\W\d\d\s\W\s([\d\s]+)\W\s\s\s\s[\s\W]'
    
    # Create a list in which there are the best clusters
    
    best_cluster = []
    with open(file,"r") as f:
        content = f.read()
        for i in re.finditer(pattern_1,content):
            best_cluster.append(i.group(1))
    
    # Cast the numbers into an integer
    
    best_cluster_int = []
    for i in best_cluster:
         best_cluster_int.append(int(i))
    
    pattern_tot = r"(?s)^(\w{3}:)[^\S\n]*(\d+)[^\S\n]*/[^\S\n]*(\d+)\n(.*?\nDOCKED: ENDMDL)"
    poses_raw = (re.findall(pattern_tot,content, re.MULTILINE))
    
    pattern_lig = r"(?s)^(\w{3}:)[^\S\n]*(\d+)[^\S\n]*/[^\S\n]*(\d+)\n(.*?\nDOCKED: TORSDOF)"
    poses_lig = (re.findall(pattern_lig,content, re.MULTILINE))
    
    def clean(pose):
        poses = []
        poses_clean = []
        
        for i in pose:
            poses.append(i[3] + "\n")
        
        for i in poses:
            poses_clean.append(re.sub("DOCKED:\s","",i))
        
        return poses_clean

    return clean(poses_lig), best_cluster_int, basename, clean(poses_raw)

def grab_table(file_pdbqt):
    '''
    
    Grab proprieties for each run of a molecule and return a the value
    
    '''
    
    with open(file_pdbqt,"r") as f:
        data = f.read()
    
    pattern_free_energy = r"(Estimated Free Energy of Binding).*?([-][\d\W]*)"
    pattern_inter = r"(Final Intermolecular Energy).*([-+][\d\W]*)"
    pattern_vdw_hbond_desolv = r"(vdW \+ Hbond \+ desolv Energy).*([-+][\d\W]*)"
    pattern_electrostatic = r"(Electrostatic Energy).*([-+][\d\W]*)"
    pattern_internal_energy = r"(Final Total Internal Energy).*([-+][\d\W]*)"
    pattern_torsional = r"(Torsional Free Energy).*([-+][\d\W]*)"
    
    def value(pattern,data):
        
        name = []
        value = []
        
        for x in re.finditer(pattern, data):
            name.append(x.group(1))
            value.append(x.group(2))
            
        value = [float(x) for x in value]
        value = tuple(value)
        table = list(zip(name,value))
        
        return value
    
    free_energy_values = value(pattern_free_energy,data)
    intermol_values = value(pattern_inter,data)
    vdw_hbond_values = value(pattern_vdw_hbond_desolv,data)
    electro_values = value(pattern_electrostatic,data)
    internal_values = value(pattern_internal_energy,data)
    torsional_values = value(pattern_torsional,data)
    
    return free_energy_values, intermol_values, vdw_hbond_values, electro_values, internal_values, torsional_values