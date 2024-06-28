import warnings
warnings.filterwarnings('ignore', 'Optimal rotation is not uniquely or poorly defined for the given sets of vectors.')

from tqdm import tqdm
import numpy as np
from pandas import DataFrame
from scipy.spatial.transform import Rotation as R

from .PDB_dotarray_volume import count_points
from .utils import timer, rnd
from .typing import *


def determine_inner(Coordinates:Points, a_pair:Tuple[str, str], a_pair_xyz:Dict[str, Point], a_pair_info:Dict[str, Contact], polar:bool) -> Surfaces:
    surface_pair = {}
    atom1 = a_pair[0]
    atom2 = a_pair[1]
    xyz: Dict[str, Point] = {}
    
    # Place atom1 at the origin
    for atom in a_pair_xyz:
        xyz[atom] = a_pair_xyz[atom] - a_pair_xyz[atom1]
    
    # Vector from atom1 to atom2
    vector: Point = xyz[atom2]
    # Calculate the rotation between this vector and the X-axis
    rotation, _ = R.align_vectors([[1, 0, 0]], [vector])
    # Rotate all points, so atom2 will fall on the X-axis, possibly in the positive direction or in the negative direction
    centers_array = np.asarray(list(xyz.values()))
    rotated_centers_array = rotation.apply(centers_array)
    # Coordinates of the points after rotation
    centers: Dict[str, Point] = {key: value for key, value in zip(xyz.keys(), rotated_centers_array)}
    
    R1: float = a_pair_info[atom1][0]
    R2: float = a_pair_info[atom2][0]
    S: DTYPE = centers[atom2][0] # The distance between the two atoms (but it may be a negative value, depending on the direction of rotation)

    # The volume represented by each point
    a_surf1 = a_pair_info[atom1][1] / len(Coordinates)
    a_surf2 = a_pair_info[atom2][1] / len(Coordinates)

    CELE = ['CA', 'C', 'CB', 'CG', 'CD', 'CE', 'CZ', 'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'] # 碳元素
    NOPS = ['N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NZ', 'NH1', 'NH2', 'O', 'OG', 'OG1', 'OD1', 'OD2', 'OE1', 'OE2', 'OH', 'P', 'SG', 'SD'] # 氮氧磷硫元素
    ELE  = {'N': ['N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NZ', 'NH1', 'NH2'], 
            'O': ['O', 'OG', 'OG1', 'OD1', 'OD2', 'OE1', 'OE2', 'OH'], 
            'P': ['P'], 
            'S': ['S', 'SG', 'SD'], 
           }
    
    # x3 and x4 are the central coordinates of two hypothetical water molecules. x4 is not necessarily used and is only assumed to exist, different from x3, when atom2 is a hydrophilic atom and the embedding of x3 into atom1 exceeds the threshold
    x4 = False
    if polar:
        atom1_type = atom1.split('+')[1]
        atom2_type = atom2.split('+')[1] # Select different embedding depths by determining whether atom2 is a hydrophilic atom.

        if atom2_type in NOPS:
            R3 = 1.4
            R4 = R2
            if S > 0:
                x3 = (S + R1 - R2) / 2 # Coordinates of the mediating water molecule
                depth  = R1 - (x3 - 1.4) # Embedding depths
                
                # The maximum embedding depths for N, O, P, and S are 0.1, 0.2, 0.3, and 0.5
                if   atom1_type in ELE['N'] and depth > 0.1: x3 = R1 + 1.4 - 0.1
                elif atom1_type in ELE['O'] and depth > 0.2: x3 = R1 + 1.4 - 0.2
                elif atom1_type in ELE['P'] and depth > 0.3: x3 = R1 + 1.4 - 0.3
                elif atom1_type in ELE['S'] and depth > 0.5: x3 = R1 + 1.4 - 0.5
                    
                depth2 = x3 + 1.4 - (S - R2)
                x4 = S - R2 - 1.4 + 0.1 if (atom2_type in ELE['N'] and depth2 > 0.1) else x3
                x4 = S - R2 - 1.4 + 0.2 if (atom2_type in ELE['O'] and depth2 > 0.2) else x3
                x4 = S - R2 - 1.4 + 0.3 if (atom2_type in ELE['P'] and depth2 > 0.3) else x3
                x4 = S - R2 - 1.4 + 0.5 if (atom2_type in ELE['S'] and depth2 > 0.5) else x3
                
            elif S < 0:
                x3 = (S - R1 + R2) / 2
                depth = R1 - (-x3 - 1.4)
                
                if   atom1_type in ELE['N'] and depth > 0.1: x3 = -R1 - 1.4 + 0.1
                elif atom1_type in ELE['O'] and depth > 0.2: x3 = -R1 - 1.4 + 0.2
                elif atom1_type in ELE['P'] and depth > 0.3: x3 = -R1 - 1.4 + 0.3
                elif atom1_type in ELE['S'] and depth > 0.5: x3 = -R1 - 1.4 + 0.5
                
                depth2 = S + R2 - (x3 - 1.4)
                x4 = S + R2 + 1.4 - 0.1 if (atom2_type in ELE['N'] and depth2 > 0.1) else x3
                x4 = S + R2 + 1.4 - 0.2 if (atom2_type in ELE['O'] and depth2 > 0.2) else x3
                x4 = S + R2 + 1.4 - 0.3 if (atom2_type in ELE['P'] and depth2 > 0.3) else x3
                x4 = S + R2 + 1.4 - 0.5 if (atom2_type in ELE['S'] and depth2 > 0.5) else x3

        else:
            if S > 0:
                x3 = (S - R2 - 1.4)
                depth = R1 - (x3 - 1.4)
                if   atom1_type in ELE['N'] and depth > 0.1: x3 = R1 + 1.4 - 0.1
                elif atom1_type in ELE['O'] and depth > 0.2: x3 = R1 + 1.4 - 0.2
                elif atom1_type in ELE['P'] and depth > 0.3: x3 = R1 + 1.4 - 0.3
                elif atom1_type in ELE['S'] and depth > 0.5: x3 = R1 + 1.4 - 0.5
                        
            elif S < 0:
                x3 = (S + R2 + 1.4)
                depth = R1 - (-x3 - 1.4)
                if   atom1_type in ELE['N'] and depth > 0.1: x3 = -R1 - 1.4 + 0.1
                elif atom1_type in ELE['O'] and depth > 0.2: x3 = -R1 - 1.4 + 0.2
                elif atom1_type in ELE['P'] and depth > 0.3: x3 = -R1 - 1.4 + 0.3
                elif atom1_type in ELE['S'] and depth > 0.5: x3 = -R1 - 1.4 + 0.5
        
        # Define the mediating water molecule as atom2
        D = S
        if atom2_type not in CELE:
            R2 = 1.4
            S  = x3       

    # Solve the equation to get the coordinates of the intersection of the two spheres on the X-axis as the threshold
    t: DTYPE = (R1**2 - R2**2 + S**2) / (2*S)
    
    atom1_coor = R1 * Coordinates
    atom2_coor = R2 * Coordinates

    # Add the coordinates of the atom center, but only the X-axis of atom2 needs to be added, because atom1 is at the origin
    atom2_coor[:, 0] += S
    
    if S > 0:   # The two atoms in a pair will have different radii and surface areas
        atom1_contacts = atom1_coor[atom1_coor[:, 0] > t] # The dot array where atom1 and atom2 contact
        atom2_contacts = atom2_coor[atom2_coor[:, 0] < t] # The dot array where atom2 and atom1 contact
    elif S < 0: # If atom2 is on the negative direction of the X-axis, then judge in reverse
        atom1_contacts = atom1_coor[atom1_coor[:, 0] < t]
        atom2_contacts = atom2_coor[atom2_coor[:, 0] > t]
    
    # Then go to calculate the number of repetitions. In this step, its own coordinates are removed
    atom1_centers = {(c, a_pair_info[c][0]): centers[c] for c in centers}; del atom1_centers[(atom1, a_pair_info[atom1][0])]
    atom2_centers = {(c, a_pair_info[c][0]): centers[c] for c in centers}; del atom2_centers[(atom2, a_pair_info[atom2][0])]
    if polar:
        atom1_centers = {(atom2, 1.4): np.asarray([S, 0, 0])}
    atom1_counts = count_points(atom1_contacts, atom1_centers)
    atom2_counts = count_points(atom2_contacts, atom2_centers)
    
    non_zero_counts1 = atom1_counts[atom1_counts != 0]
    non_zero_counts2 = atom2_counts[atom2_counts != 0]
    surf1: DTYPE = (1 / non_zero_counts1).sum()   # Add atom2 itself
    surf2: DTYPE = (1 / non_zero_counts2).sum()
    
    # Put the distance in
    dist   = rnd(np.abs(S))
    f_surf = rnd(surf1 * a_surf1)
    r_surf = rnd(surf2 * a_surf2)
    type1  = a_pair_info[atom1][-1]
    type2  = a_pair_info[atom2][-1]

    if x4:
        t2: DTYPE = (R4**2 - D**2 + x4**2 - 1.4**2) / (2*(x4 - D))
        atom3_coor = R3 * Coordinates
        atom4_coor = R4 * Coordinates
        atom3_coor[:, 0] += x4
        atom4_coor[:, 0] += D
        if D > 0:
            atom3_contacts = atom3_coor[atom3_coor[:, 0] > t2]
            atom4_contacts = atom4_coor[atom4_coor[:, 0] < t2]
        elif D < 0:
            atom3_contacts = atom3_coor[atom3_coor[:, 0] < t2]
            atom4_contacts = atom4_coor[atom4_coor[:, 0] > t2]
        atom4_centers = {(atom1, 1.4): np.asarray([x4, 0, 0])}
        atom4_counts = count_points(atom4_contacts, atom4_centers)
        non_zero_counts4 = atom4_counts[atom4_counts != 0]
        surf4: DTYPE = (1 / non_zero_counts4).sum()
        r_surf = rnd(surf4 * a_surf2)
        f_surf = (f_surf + r_surf) / 2

    if polar:
        dist = rnd(np.abs(D))
        surface_pair[a_pair]       = (f_surf, dist, type1, type2)
    else:
        surface_pair[a_pair]       = (f_surf, dist, type1, type2)
        surface_pair[a_pair[::-1]] = (r_surf, dist, type2, type1)
    
    return surface_pair


def pdb_dotarray_surface(ref_fp:Path, contact_df:DataFrame, atom_pairs:dict, polar:bool, disable_print=False) -> Surfaces:
    @timer(disable_print=disable_print)
    def count_surface():
        # Coordinates of the atom center
        atoms = contact_df['Name']
        atom_coords = np.asarray([
            contact_df['x'].array, 
            contact_df['y'].array,
            contact_df['z'].array,
        ], dtype=DTYPE).T
        contacts_center: Dict[str, Point] = {i: j for i, j in zip(atoms, atom_coords)}
        contacts_dict: Dict[str, Contact] = contact_df[['Name', 'R', 'Surf', 'Type']].set_index('Name').T.to_dict('list')

        # Import the coordinates of the points
        Coordinates: Points = np.loadtxt(ref_fp, dtype=DTYPE)
        Pairs_Surface = []
        for a_pair in tqdm(atom_pairs, disable=disable_print):
            # A pair of atoms, and all atoms that contact them
            a_atom_pairs = atom_pairs[a_pair]
            a_pair_xyz  = {name: contacts_center[name] for name in a_atom_pairs if name in contacts_center}     # Coordinates
            a_pair_info = {name: contacts_dict  [name] for name in a_atom_pairs if name in contacts_dict  }     # Radius, surface, type
            
            surface_pair = determine_inner(Coordinates, a_pair, a_pair_xyz, a_pair_info, polar)
            Pairs_Surface.append(surface_pair)
        
        surface = {} # Use update() to quickly merge dictionaries
        for pair in Pairs_Surface:
            surface.update(pair)  

        return surface
    
    return count_surface()
