The script requires:
- pandas = conda install -c anaconda pandas
- rdkit = conda install -c conda-forge rdkit
- openbabel = conda install -c conda-forge openbabel

There are two scripts:

- table_and_rdkit = it will do a conformational search on each ligand, than it will calculate the internal conformational energy of each pose and make a delta (ligand energy - global minimum). Also it will grab the autodock parameters and add everything to a tidy excel file (if you do prefer a csv, you can just edit the last lines).
As input are needed the output files from autodock GPU (.dlg) and the smiles of the ligand (both in the working directory), be careful that each smile has to be in a different file since the script will grab the filename. Also the .smiles and the .dlg files have to have the same name

- table_autodock = it will just write a table with the information provided by autodock for each ligand. As input are just needed the files from autodock GPU (.dlg)

The chosen poses are the best poses of each cluster, clusterized by the parameters of autodock.
