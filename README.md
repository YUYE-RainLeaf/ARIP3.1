# ARIP3

A Python software for quantitative calculation of residue interactions in proteins or nucleic acids, supporting platforms such as Windows/Linux/MacOS.

It can simultaneously calculate contact area and volume, and supports analysis of PDB files containing multiple MODELs.
Supports most non-standard residues and ligands.

----

### Getting Started

- Ensure your Python version is 3.6 or higher. You can check your Python version by running `python --version`
- Download and use immediately
  - `pip install -r requirements.txt`
  - `python3 run.py your/file/or/dir`
  
- Supports input of files or paths
  - The input file must be in PDB format, both compressed and uncompressed are acceptable, with or without suffix
  - The input path should only contain correctly formatted files and no secondary paths
  - If a path is input, multithreading will be used automatically
  - Water molecules will be automatically discarded when reading non-standard residues

### Optional Parameters

- Custom output path
  - `python3 run.py your/file/or/dir -o your/output/path`
  
- Enhanced precision mode, more time-consuming
  - `python3 run.py your/file/or/dir -e`
  
- Use the lower cutoff, contacts with area and volume less than the cutoff will be discarded. The default is 0.5 and 0.2
  - `python3 run.py your/file/or/dir -c` # Contacts with area less than 0.5 and volume less than 0.2 will be discarded
  - `python3 run.py your/file/or/dir -c 1.0 0.5` # Contacts with area less than 1.0 and volume less than 0.5 will be discarded
  - Note, this will apply to all proteins, nucleic acids, non-standard residues, and ligands

- Use a custom distance to determine whether two atoms are in contact. Contact is determined when the minimum distance between two atoms is less than the custom distance. The default is 0
- `python3 run.py your/file/or/dir -d` # Contact is determined when the minimum distance between two atoms is less than 0
- `python3 run.py your/file/or/dir -d 0.5` # Contact is determined when the minimum distance between two atoms is less than 0.5
- The minimum distance refers to the smallest distance between two atoms when they are considered as spheres with their vdW

- Custom multithreading quantity
  - `python3 run.py your/file/or/dir -t threads_number`
  - If the -t parameter is not specified, or the number of threads entered is greater than the number of CPUs, multithreading will be set automatically according to the number of CPUs
  
- Calculate area and volume according to the algorithm of hydrophilic atoms interactions mediated by water molecule
  - `python3 run.py your/file/or/dir -p`
  - This mode only calculates atom pairs where at least one atom is N, O, P, or S, with the nearest distance between two atoms > 0 and < 2.8
  - The maximum embedding depths for N, O, P, and S are 0.1, 0.2, 0.3, and 0.5, respectively. If exceeded, it will automatically be corrected to these thresholds
  - If only one of the two atoms is N, O, P, or S, the calculated area and volume are the contact area and volume of that atom with one water molecule
  - If both atoms are N, O, P, or S, the calculated area and volume are the average of the contact areas and volumes of these two atoms with one water molecule each
  - Since the values are usually small, this mode further increases the number of points in the volume calculation grid to achieve higher accuracy
  - In this mode, custom distance criteria for determining contact between two atoms cannot be used, and the effect of multiple atoms in contact at the same time will not be considered
  - When two carbon atoms come into contact, only direct contact is considered as contact, without the mediation of water molecules

- Calculate SASA (Solvent Accessible Surface Area)
  - `python3 run.py your/file/or/dir -a`
  - This mode calculates the surface area on each atom within a residue that is not in contact with any other atom, and then sums them up
  - The results of the SASA calculation are in the CSV file with prefix '_SUM'.
  
- Calculate volume based on atomic overlap weighted algorithm
  - `python3 run.py your/file/or/dir -w`
  - The more overlapping atoms, the greater the volume weight
  - Specifically, the volume represented by each point is multiplied by the number of atoms that contain this point, and if it is located in N atoms at the same time, it will be ×N
  
- Save the output results in .gz compressed format to save storage space
  - `python3 run.py your/file/or/dir -z`
  - The created folder name will have a '_z' suffix to distinguish it

- Output .csv file for each residue
  - `python3 run.py your/file/or/dir -r`
  - This may consume storage space and reduce write speed during large-scale runs

- Calculate only surface area and other parameters, excluding volume
  - `python3 run.py your/file/or/dir -s`
  - This can be used when volume values are not required, significantly improving speed 

### Output Style

##### Folder
- For each successfully analyzed PDB file, a new folder with the same name will be created to store the analysis results
- If the PDB file contains multiple MODELs, subfolders will be created under the same name folder to store them separately

##### Interaction Types

- Protein
  - Non-covalent interaction (NC, Non-Covalent)
	- HB: Hydrogen Bond
    - AROM: Aromatic-Aromatic Contact
    - PHOB: Hydrophobic-Hydrophobic Contact
	- DC: Distabilizing Contact
	- OTHER: Other van der Waals interactions
  - Covalent interaction (Cova, Covalent)
    - SS: Disulfide Bond
	- PB: Peptide Bond
	
- Nucleic Acid
  - Non-covalent interaction
    - DD: DNA-DNA Contact
	- DR: DNA-RNA Contact
	- RR: RNA-RNA Contact
  - Covalent interaction
    - PD: Phosphodiester Bond
	
- Other
  - SASA: Solvent Accessible Surface Area, surface area on the residue not in contact with any atom
  - Surf: Surface area, the sum of the contact surface area of each atom in the residue
  - Volu: Volume, contact volume
  - AOWV: Atomic Overlap Weighted Volume
  - UNDEF: Undefined, interactions related to non-standard residues that have not been judged

##### Summary Information
- Files with preffix '_ALL' contain all atom-atom contact information, including atom type, contact distance, surface, volume, etc.
- Files with preffix '_RES' contain summarized information on pairwise contacts between residues.
- Files with preffix '_SUM' contain each residue's dihedral angles (proteins only), SASA (if available), covalent and non-covalent contact areas, and contact volumes, etc.

##### Success Information
- If the PDB file only has one MODEL
  - "The PDB {name} run OK, time cost: {time}s"

- If the PDB file contains multiple MODELs
  - "The PDB {name}MODEL{num} run OK, time cost: {time}s"

#### Error Information

- The input is not a valid file or path
  - "{input} is not a valid file or directory"

- The cutoff is specified, but no correctly formatted numbers are entered
  - "Lower cutoff must be TWO values like 1.0 0.5, or leave it blank to use the default values 0.5 0.2"

- The PDB file is incomplete
  - "The file {name} may be corrupted or contain incomplete MODEL"

- The PDB format is incorrect
  - "The file {name} is unsupported format"

- Insufficient memory or the input file does not contain valid atoms
  - "Required memory: {required_memory_GB} GB. Skipping file {name} due to insufficient memory or no valid atoms"

- Unable to analyze for other reasons
  - "The file {name} cannot be analyzed, perhaps it contains unsupported format, or no valid atoms"
  
### Visualization

- Provides a simple function to visualize the spatial position of heavy atoms
  - `python3 vis.py your/file`
  
  - Supports input of PDB files or XYZ files
    - PDB files, both compressed and uncompressed are acceptable, with or without suffix
    - Does not support input of paths
    - For PDB containing multiple MODELs, pictures of each MODEL can be output
    - If an unsupported format is input, it will output `unsupported file: {name}`

- For structures that have been analyzed, heatmaps of residue contacts can be drawn for PDBs with multiple models
  - `python3 heatmap.py your/dir column`
  - The input path should contain the ARIP analysis results
    - Only applicable to PDBs with multiple models
    - The input column is the column you wish to plot in the heatmap, generally one or more of the following:
      - Surface_max, Surface_min, Surface_range, Surface_mean, Surface_rmsd, Surface_cv, Volume_max, Volume_min, Volume_range, Volume_mean, Volume_rmsd, Volume_cv, AOWV_max, AOWV_min, AOWV_range, AOWV_mean, AOWV_rmsd, AOWV_cv
      - If no input column, a heatmap will be drawn for each column by default
    - The x-axis and y-axis of the heatmap are arranged in the order of the residues, starting from A1.
    - If an unsupported format is input, it will output `The heatmap could not be successfully generated due to an error`  
  
### Contact me
Email: xiangtao312@outlook.com

----
2024/06/28