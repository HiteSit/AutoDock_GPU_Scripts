#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob, os, shutil
import re
import pdb

import pandas as pd
import numpy as np

from openbabel import pybel, openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import chain


# In[2]:


from functions import parse_mol as parse
from functions import conformer_gen as conf


# In[ ]:


# Choose work directory
work_dir = input("Choose your workdir: ")


# In[ ]:


# Define choice-functions
def choice_num_conf():
    choice = "wrong"
    
    while choice.isdigit() == False:
        choice = input("\nChoose number of conformations to generate: ")
        
        if choice.isdigit() == False:
            print("\nYou must provide a number")
    
    return int(choice)

def choice_max_iter():
    choice = "wrong"
    
    while choice.isdigit() == False:
        choice = input("\nChoose max interactions: ")
        
        if choice.isdigit() == False:
            print("\nYou must provide a number")
    
    return int(choice)


# In[ ]:


# Move to the working directory
os.chdir(work_dir)


# In[ ]:


# Write files full name & files basename
files = glob.glob("*.dlg")
basenames = []
for file in files:
    basenames.append(file.split(".")[0])


# In[8]:


# Grab the smiles of the ligands
smile_files = glob.glob("*.smiles")

# Generate some conformer for each ligand, minimize all the conformations and write the minimum in a list of tuples
print("\nGenerating the conformation for each ligand")
num_conf = choice_num_conf()
max_iter = choice_max_iter()


smiles = []
min_confs = []
for file in smile_files:
    with open(file, "r") as f:
        data = f.readline()
        smiles.append(data)
    
    min_confs.append(conf.gen_conf(data, num_conf, max_iter))


# In[ ]:


tup_1 = list(zip(files,basenames))
conv_for_pose = []

print("Writing pdbqt files")

for file, basename in tup_1:
    
    # Grab string 100 poses (list len=100)
    data = parse.grab_mol_clean(file)[0]
    data_save = parse.grab_mol_clean(file)[3]
    
    # Grab clusters run
    clusters = parse.grab_mol_clean(file)[1]
    
    # Convert only the best poses for each run a append them to a list ligand-indipendent
    for num in clusters:
        conv_for_pose.append(conf.convert_to_mol2(data)[num-1])
    
    # Write PDBQT
    with open(basename + ".pdbqt", "w") as f:
        for num in clusters:
            f.write(data_save[num-1])
    
    pdbqt = glob.glob("*.pdbqt")


# In[ ]:


# Calculate conformational energy for each pose in the previusly made list (still ligand-indipendent)
energy_for_run = []

for pose in conv_for_pose:
    energy_for_run.append(conf.calc_conf_energy(pose))


# In[ ]:


# Creare a tuple with all the run number in a ligand-indipendent way
runs_nest = []

for file in files:
    runs_nest.append(parse.grab_mol_clean(file)[1])

runs = tuple(chain.from_iterable(runs_nest))


# In[ ]:


# Link the molecule with the global minimum conformer
min_conf_list = []

for f,l in min_confs:
    l_round = round(l,2)
    min_conf_list.append(l_round)

mol_and_min_conf = list(zip(files, min_conf_list))


# In[ ]:


# Link the molecule and the conf_energy with the number of runs so that every list has the same lenght
mol_for_run = []
conf_for_run = []

for file, conf in mol_and_min_conf:
    temp = parse.grab_mol_clean(file)[1]
    
    for t in temp:
        mol_for_run.append(file)
        conf_for_run.append(conf)


# In[ ]:


# Zip everything
mol_conf_for_run = list(zip(mol_for_run,conf_for_run,energy_for_run,runs))


# In[ ]:


# Calculate Delta conformational energy
final_list = []

for a,b,c,d in mol_conf_for_run:
    
    temp = []
    
    temp.append(a)
    temp.append(round(c-b,2))
    temp.append(d)
    
    temp_tuple = tuple(temp)
    final_list.append(temp_tuple)
    
final_tuple = tuple(final_list)


# In[ ]:


# Make index data frame
index = pd.MultiIndex.from_tuples(final_tuple, names=["Molecule","Delta_Conf_Energy","Run"])
df = pd.DataFrame(index=index)


# In[ ]:


# Prepare columns
energy = []
intermol_energy = []
vdw_hbond_energy = []
electrostatic_energy = []
internal_energy = []
torsional_energy = []

for file in pdbqt:
    energy.append(parse.grab_table(file)[0])
    intermol_energy.append(parse.grab_table(file)[1])
    vdw_hbond_energy.append(parse.grab_table(file)[2])
    electrostatic_energy.append(parse.grab_table(file)[3])
    internal_energy.append(parse.grab_table(file)[4])
    torsional_energy.append(parse.grab_table(file)[5])

energy = tuple(chain.from_iterable(energy))
intermol_energy = tuple(chain.from_iterable(intermol_energy))
vdw_hbond_energy = tuple(chain.from_iterable(vdw_hbond_energy))
electrostatic_energy = tuple(chain.from_iterable(electrostatic_energy))
internal_energy = tuple(chain.from_iterable(internal_energy))
torsional_energy = tuple(chain.from_iterable(torsional_energy))


# In[ ]:


# Append to index data frame
df_final = df.assign(
    Free_Energy_Of_Binding = energy,
    Intermolecular_Energy = intermol_energy,
    VdW_HBond_Desolv_Energy = vdw_hbond_energy,
    Electrostatic_Energy = electrostatic_energy,
    Internal_Energy = internal_energy,
    Torsional_Free_Energy = torsional_energy
)


# In[ ]:


# Create a new dataframe sorting the value by dock score
df_sort = df_final.sort_values("Free_Energy_Of_Binding",ascending=True)


# In[ ]:


# Order the files in the right directories
try:
    os.mkdir("pdbqt")
except:
    pass

for file in pdbqt:
    dest = "pdbqt/" + file
    shutil.move(file, dest)

os.chdir("pdbqt/")


# In[ ]:


# Write excel files
df_sort.to_excel("Pose_Table_Sort.xlsx")
df_final.to_excel("Pose_Table.xlsx")

