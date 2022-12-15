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


# In[4]:


# Choose work directory
work_dir = input("Choose your workdir: ")

# Move to the working directory
os.chdir(work_dir)


# In[5]:


# Write files full name & files basename
files = glob.glob("*.dlg")
basenames = []
for file in files:
    basenames.append(file.split(".")[0])


# In[6]:


tup_1 = list(zip(files,basenames))
conv_for_pose = []
for file, basename in tup_1:
    
    # Grab string 100 poses (list len=100)
    data_save = parse.grab_mol_clean(file)[3]
    
    # Grab clusters run
    clusters = parse.grab_mol_clean(file)[1]
    
    # Write PDBQT
    with open(basename + ".pdbqt", "w") as f:
        for num in clusters:
            f.write(data_save[num-1])
    
    pdbqt = glob.glob("*.pdbqt")


# In[7]:


# Creare a tuple with all the run number in a ligand-indipendent way
runs_nest = []

for file in files:
    runs_nest.append(parse.grab_mol_clean(file)[1])

runs = tuple(chain.from_iterable(runs_nest))


# In[8]:


mol_for_run = []

for file in files:
    temp = parse.grab_mol_clean(file)[1]
    
    for t in temp:
        mol_for_run.append(file)


# In[9]:


# Zip everything
mol_conf_for_run = list(zip(mol_for_run,runs))


# In[10]:


# Make index data frame
index = pd.MultiIndex.from_tuples(mol_conf_for_run, names=["Molecule","Run"])
df = pd.DataFrame(index=index)


# In[11]:


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


# In[12]:


# Append to index data frame
df_final = df.assign(
    Free_Energy_Of_Binding = energy,
    Intermolecular_Energy = intermol_energy,
    VdW_HBond_Desolv_Energy = vdw_hbond_energy,
    Electrostatic_Energy = electrostatic_energy,
    Internal_Energy = internal_energy,
    Torsional_Free_Energy = torsional_energy
)


# In[13]:


# Create a new dataframe sorting the value by dock score
df_sort = df_final.sort_values("Free_Energy_Of_Binding",ascending=True)


# In[14]:


# Order the files in the right directories
try:
    os.mkdir("pdbqt")
except:
    pass

for file in pdbqt:
    dest = "pdbqt/" + file
    shutil.move(file, dest)

os.chdir("pdbqt/")


# In[15]:


# Write excel files
df_sort.to_excel("Pose_Table_Sort.xlsx")
df_final.to_excel("Pose_Table.xlsx")

