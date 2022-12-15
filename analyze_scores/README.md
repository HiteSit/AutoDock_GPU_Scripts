The program requires:

- pandas = conda install -c anaconda pandas
- rdkit = conda install -c conda-forge rdkit
- openbabel = conda install -c conda-forge openbabel

There are two scripts:

- table_and_rdkit = it will do a conformational search on each ligand, than it will calculate the internal conformational energy of each pose and make a delta (ligand energy - global minimum). Also it will grab the autodock parameters and add everything to a tidy excel file (if you do prefer a csv, you can just edit the last lines). This script will ask for the number of conformer to generate, the maximum interaction to consider and the working directory
  - As input are needed the output files from autodock GPU (.dlg) and the smiles of the ligand (both in the working directory), be careful that each smile has to be in a different file since the script will grab the filename. Also the .smiles and the .dlg files have to have the same name.

- table_autodock = it will just write a table (excel) with the informations provided by autodock for each ligand. This script will just ask for the working directory.
  - As input are just needed the files from autodock GPU (.dlg).

The poses selected by the script are the ones compliance with the table in the end of each autodock files. So the each pose is the best pose of its cluster.
