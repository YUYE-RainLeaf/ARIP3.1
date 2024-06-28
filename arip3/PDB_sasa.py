import warnings
warnings.filterwarnings('ignore', 'Optimal rotation is not uniquely or poorly defined for the given sets of vectors.')

from tqdm import tqdm
import copy
import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.api.types import is_numeric_dtype
from scipy.spatial.transform import Rotation as R

from .PDB_constants import *
from .PDB_dotarray_volume import count_points
from .utils import timer, rnd
from .typing import *


def determine_inner(Coordinates:np.ndarray, R1:float, a_sasa, a_atom_list:List[str], a_atom_xyz:Dict[str, Point], a_atom_info:Dict[str, Contact]) -> float:
    atom_coor = R1 * Coordinates
    for atom in a_atom_list:
        R2 = a_atom_info[atom][0]
        # Vector from atom1 to atom2
        vector: Point = a_atom_xyz[atom]
        # Calculate the rotation between this vector and the X-axis
        rotation, _ = R.align_vectors([[1, 0, 0]], [vector])
        # Rotate the coordinates in atom_coor and a_atom_xyz
        atom_coor = rotation.apply(atom_coor)
        for atom_key in a_atom_xyz:
            a_atom_xyz[atom_key] = rotation.apply(a_atom_xyz[atom_key][np.newaxis, :])[0]
        
        S: DTYPE = a_atom_xyz[atom][0] # The distance between the two atoms (but it may be a negative value, depending on the direction of rotation)
        t: DTYPE = (R1**2 - R2**2 + S**2) / (2*S) # Solve the equation to get the coordinates of the intersection of the two spheres on the X-axis as the threshold
        
        # Retain points where the x-coordinate is less than t
        if   S > 0:
            atom_coor = atom_coor[atom_coor[:, 0] < t]
        elif S < 0:
            atom_coor = atom_coor[atom_coor[:, 0] > t]

    sasa = a_sasa * len(atom_coor)
    
    return sasa


def pdb_dotarray_sasa(ref_fp:Path, atom_model:AtomModel, disable_print=False) -> SASA:
    def parse_atom_model() -> DataFrame:
        # Create a dictionary to convert to DataFrame
        atom_info = {
            'Chain':   [],      # str
            'ResName': [],      # str
            'Atom':    [],      # str
            'x':       [],      # float
            'y':       [],      # float
            'z':       [],      # float
            'R':       [],      # float
            'Surf':    [],      # float
        }

        for atom in atom_model:
            Res = atom[17:20].lstrip()  # Residue name
            Num = str(int(atom[22:26])) # Residue number
            Ana = atom[12:16].strip()   # Atom name
            Ele = atom[76:78].lstrip()  # Element type
            
            # Residue, atom, coordinates, radius, type, surface, volume
            atom_info['Chain']  .append(atom[21] + Num)
            atom_info['ResName'].append(Res)
            atom_info['Atom']   .append(Ana)
            atom_info['x']      .append(float(atom[30:38]))
            atom_info['y']      .append(float(atom[38:46]))
            atom_info['z']      .append(float(atom[46:54]))
            
            if Res in Radius: # Standard residue
                if Res in AA and Ana in Radius[Res]: # Amino acid
                    __ = Radius[Res][Ana]
                    atom_info['R']   .append(__[0])
                    atom_info['Surf'].append(__[2])
                        
                if Res in NT and Ana[0] in Radius[Res]: # Nucleotide
                    __ = Radius[Res][Ele]
                    atom_info['R']   .append(__[0])
                    atom_info['Surf'].append(__[2])
            
            elif Res in Radius['UNDEF']: # Non-standard residue
                __ = Radius['UNDEF'][Ele]
                atom_info['R']   .append(__[0])
                atom_info['Surf'].append(__[2])
            else:
                __ = Radius['UNDEF']['X']
                atom_info['R']   .append(__[0])
                atom_info['Surf'].append(__[2])
        
        atom_df = pd.DataFrame.from_dict(atom_info)
        
        # NOTE: Add the radius of the water molecule
        atom_df['R'] += Radius_H2O
        
        for col in atom_df.columns:
            if is_numeric_dtype(atom_df[col]): # Convert all values to np.float64 type
                atom_df[col] = atom_df[col].astype(DTYPE, copy=False)
        
        # Extract rows in the ResName column in the AA or NA dictionary
        aa = atom_df[atom_df['ResName'].isin(AA)]
        nt = atom_df[atom_df['ResName'].isin(NT)]
        ns = atom_df[~atom_df['ResName'].isin(AA) & ~atom_df['ResName'].isin(NT)]
        
        # Merge Chain and ResName columns
        aa_tmp = aa.copy()
        nt_tmp = nt.copy()
        ns_tmp = ns.copy()
        aa_tmp['ResName'] = aa_tmp['ResName'].map(AA)
        nt_tmp['ResName'] = nt_tmp['ResName'].map(NT)
        
        # Use a different connector for easy replacement
        aa_tmp['Name'] = aa_tmp['Chain'] + ',' + aa_tmp['ResName'] + '+' + aa_tmp['Atom']
        nt_tmp['Name'] = nt_tmp['Chain'] + ',' + nt_tmp['ResName'] + '+' + nt_tmp['Atom']
        ns_tmp['Name'] = ns_tmp['Chain'] + ';' + ns_tmp['ResName'] + '+' + ns_tmp['Atom']
        
        aa_df = aa_tmp[['Name', 'x', 'y', 'z', 'R', 'Surf']]
        nt_df = nt_tmp[['Name', 'x', 'y', 'z', 'R', 'Surf']]
        ns_df = ns_tmp[['Name', 'x', 'y', 'z', 'R', 'Surf']]
            
        atom_df = pd.concat([aa_df, nt_df, ns_df], axis=0)
        
        return atom_df    
    
    @timer(disable_print=disable_print)
    def count_sasa():
        atom_df = parse_atom_model()
        # Coordinates of the atom center
        atoms = atom_df['Name']
        atom_coords = np.asarray([
            atom_df['x'].array, 
            atom_df['y'].array,
            atom_df['z'].array,
        ], dtype=DTYPE).T
        contacts_center: Dict[str, Point] = {i: j for i, j in zip(atoms, atom_coords)}
        contacts_dict: Dict[str, Contact] = atom_df[['Name', 'R', 'Surf']].set_index('Name').T.to_dict('list')

        # Build into atom pairs
        dists = np.linalg.norm(atom_coords[:, np.newaxis] - atom_coords, axis=2)
        eps = 0.00001
        radii = np.asarray(atom_df['R'])
        contact_map = {
            atom1: atoms[(dists[i] > eps) & (dists[i] < radii[i] + radii)]
                for i, (atom1, radius) in enumerate(zip(atoms, radii))
        }
        contact_dict = {} # Build into atom pairs
        for key, values in contact_map.items():
            residue_name = key.split('+')[0]
            if residue_name not in contact_dict:
                contact_dict[residue_name] = {}
            contact_dict[residue_name][key] = values.tolist()
        
        sasa = {}
        # Import the coordinates of the points
        Coordinates: Points = np.loadtxt(ref_fp, dtype=DTYPE)

        for a_residue in tqdm(contact_dict, disable=disable_print):
            # A residue, and information of all atoms that contact other atoms
            Residue_SASA = 0
            for a_atom in contact_dict[a_residue]:
                a_atom_list = contact_dict[a_residue][a_atom]
                a_atom_xyz  = {name: contacts_center[name] for name in a_atom_list if name in contacts_center}     # Coordinates
                a_atom_info = {name: contacts_dict  [name] for name in a_atom_list if name in contacts_dict  }     # Radius, surface, type
                a_atom_xyz_copy = copy.deepcopy(a_atom_xyz) # deepcopy
                
                for atom in a_atom_xyz_copy:
                    a_atom_xyz_copy[atom] -= contacts_center[a_atom]
                R1 = contacts_dict[a_atom][0]
                a_sasa = contacts_dict[a_atom][1] / len(Coordinates)
                
                sasa_atom = determine_inner(Coordinates, R1, a_sasa, a_atom_list, a_atom_xyz_copy, a_atom_info)
                Residue_SASA += sasa_atom
        
            sasa[a_residue] = Residue_SASA

        return sasa
    
    return count_sasa()
