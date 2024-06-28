import warnings
warnings.filterwarnings('ignore', 'Optimal rotation is not uniquely or poorly defined for the given sets of vectors.')

from tqdm import tqdm
import numpy as np
from scipy.spatial.transform import Rotation as R

from .utils import timer, rnd
from .typing import *


# Define a function to generate a number of points
def gen_points(x_range:Range, y_range:Range, z_range:Range, interval:float=0.2) -> Optional[Points]:
    # Calculate the size of the point array
    size_x = int((x_range[1] - x_range[0]) / interval) + 1 
    size_y = int((y_range[1] - y_range[0]) / interval) + 1
    size_z = int((z_range[1] - z_range[0]) / interval) + 1
    sizes = [size_x, size_y, size_z]
    if 1 in sizes: return

    # Build a coordinate point array: [[0, 0, 0], [0, 0, 1], ..., [nX-1, xY-1, xZ-1]]
    zeros = np.zeros(sizes)                         # [nX, xY, xZ]
    grid = np.transpose(np.where(zeros == 0.0))     # [nX*xY*xZ, D=3]
    grid = grid.astype(DTYPE)

    # Remap data: coord_rang: [0, nAxis] => v_range: (vmin, vmax)
    def map_data(idx:int, size:int, vrng:Range):
        grid[:, idx] = (grid[:, idx] / (size - 1)) * (vrng[1] - vrng[0]) + vrng[0]

    map_data(0, size_x, x_range)
    map_data(1, size_y, y_range)
    map_data(2, size_z, z_range)

    return grid


def count_points(atom_contacts:Points, atom_centers:Dict[Tuple[str, float], Point]) -> ndarray:
    # Calculate the distance between all points and the atom center
    atom_info = np.asarray(list(atom_centers.values()))                # Atom coordinates
    radii = np.asarray([radius for _, radius in atom_centers.keys()])  # Atom radius
    d = atom_info[:, np.newaxis] - atom_contacts
    dists = np.linalg.norm(d, axis=2) # This step is slow, but it seems to be the only solution

    # Find points with a distance less than or equal to 'r'
    within_r = dists <= radii[:, np.newaxis]
    
    # Calculate how many atoms each point is inside. No need for atom coordinates at all
    atom_counts = np.sum(within_r, axis=0)
    
    return atom_counts


def determine_inner(a_pair:Tuple[str, str], a_pair_xyz:Dict[str, Point], a_pair_info:Dict[str, float], polar:bool, weighted:bool, interval:float=0.2) -> Tuple[dict, dict]:
    volume_pair = {}
    xyz = {}
    atom1, atom2 = a_pair
    
    # Place atom1 at the origin
    for atom in a_pair_xyz:
        xyz[atom] = a_pair_xyz[atom] - a_pair_xyz[atom1]
    
    # Vector from atom2 to atom1
    vector = xyz[atom2]
    # Calculate the rotation between this vector and the X-axis
    rotation, _ = R.align_vectors([[1, 0, 0]], [vector])
    # Rotate all points, so atom2 will fall on the X-axis, possibly in the positive direction or in the negative direction
    centers_array = np.asarray(list(xyz.values()))
    rotated_centers_array = rotation.apply(centers_array)
    # Coordinates of the points after rotation
    centers = {key: value for key, value in zip(xyz.keys(), rotated_centers_array)}
    
    R1 = DTYPE(a_pair_info[atom1])
    R2 = DTYPE(a_pair_info[atom2])
    S  = DTYPE(centers[atom2][0])    # The distance between the two atoms (but it may be a negative value, depending on the direction of rotation)
    
    # Determine the range of the point array
    if R1 >= R2:
        y_max = + R1
        y_min = - R1
        z_max = + R1
        z_min = - R1
    elif R2 > R1:
        y_max = + R2
        y_min = - R2
        z_max = + R2
        z_min = - R2
        
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

            two_center2 = np.asarray([centers[atom2], [x4, 0, 0]])

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
        
        if S > 0:
            x_max = 0 + R1
            x_min = x3 - 1.4
            if x4:
                x_max2 = x4 + 1.4
                x_min2 = S - R2
        elif S < 0 :
            x_max = x3 + 1.4
            x_min = 0 - R1
            if x4:
                x_max2 = S + R2
                x_min2 = x4 - 1.4
                
        # Define the mediating water molecule as atom2
        if atom2_type not in CELE:
            R2 = 1.4
            
        # Coordinates of atom1 and atom2
        two_center = np.asarray([centers[atom1], [x3, 0, 0]])
        
        if atom1_type in CELE and atom2_type in CELE:
            if S > 0:
                x_max = 0 + R1
                x_min = S - R2
            elif S < 0:
                x_max = S + R2
                x_min = 0 - R1
            two_center = np.asarray([centers[atom1], centers[atom2]])
            # When two carbon atoms come into contact, only direct contact is considered as contact, without the mediation of water molecules.
    
    else:
        if S > 0:
            x_max = 0 + R1
            x_min = S - R2
        elif S < 0:
            x_max = S + R2
            x_min = 0 - R1
        two_center = np.asarray([centers[atom1], centers[atom2]])    
    
    # Calculate the volume and number of points in the point array, +0.1 is to avoid points too close to the edge being inaccurate
    eps = 0.1
    if x_max < x_min: return {}
    cube_volu: DTYPE = rnd((x_max - x_min + eps) * (y_max - y_min + eps) * (z_max - z_min + eps), 9)

    # Specify the range of each dimension
    x_range = (rnd(x_min - eps / 2), rnd(x_max + eps / 2))  
    y_range = (rnd(y_min - eps / 2), rnd(y_max + eps / 2))
    z_range = (rnd(z_min - eps / 2), rnd(z_max + eps / 2))

    # Generate point array coordinates. If there is no more than 1 layer of point array in a certain dimension, discard it directly
    points: Points = gen_points(x_range, y_range, z_range, interval=interval)
    if points is None: return {}
    
    # Distance from the point array to atom1 and atom2
    dists = np.linalg.norm(points[:, np.newaxis] - two_center, axis=2)
    
    # Points within atom1 and atom2
    in_contact = points[(dists[:, 0] <= R1) & (dists[:, 1] <= R2)]
    atom1_centers = {(c, a_pair_info[c]): centers[c] for c in centers} ; del atom1_centers[(atom1, a_pair_info[atom1])]
    if polar and not (atom1_type in CELE and atom2_type in CELE):
        atom1_centers = {(atom2, 1.4): np.asarray([x3, 0, 0])}
    atom1_counts = count_points(in_contact, atom1_centers)

    # The volume represented by each point
    a_point = cube_volu / len(points)
    
    # Finally calculate the contact volume
    # Ordinary calculation
    non_zero_counts = atom1_counts[atom1_counts != 0]
    volu: DTYPE = (1 / non_zero_counts).sum() # No need to +1, just manage the number of contacts with other atoms
    volume = rnd(volu * a_point)
    
    # If the last three digits after the decimal point are all 0, ignore it
    if volume == 0: return {}
    if x4: # The volume is the average of the contact volumes of atom1 and atom2 with the water molecule, respectively
        cube_volu2: DTYPE = rnd((x_max2 - x_min2 + eps) * (y_max - y_min + eps) * (z_max - z_min + eps), 9)
        x_range2 = (rnd(x_min2 - eps / 2), rnd(x_max2 + eps / 2))  
        points2: Points = gen_points(x_range2, y_range, z_range, interval=interval)
        if points2 is None: return {}
        dists2 = np.linalg.norm(points2[:, np.newaxis] - two_center2, axis=2)
        in_contact2 = points2[(dists2[:, 0] <= R4) & (dists2[:, 1] <= R3)]
        atom1_centers2 = {(atom1, 1.4): np.asarray([x4, 0, 0])}
        atom1_counts2 = count_points(in_contact2, atom1_centers2)
        a_point2 = cube_volu2 / len(points2)
        non_zero_counts2 = atom1_counts2[atom1_counts2 != 0]
        volu2: DTYPE = (1 / non_zero_counts2).sum()
        volume2 = rnd(volu2 * a_point2)
        volume = (volume + volume2) / 2
        
    # Atomic Overlap Weighted Algorithm
    if weighted: 
        edvolu: DTYPE = (1 + atom1_counts).sum() # +1 is because the atom itself is not counted
        edvolume = rnd(edvolu * a_point)
        
        if x4:
            edvolu2: DTYPE = (1 + non_zero_counts2).sum()
            edvolume2 = rnd(edvolu2 * a_point)
            edvolume = (edvolume + edvolume2) / 2
            
        # Add the reverse atom pair as well, for easy processing later
        volume_pair[a_pair]       = [volume, edvolume]
        if polar == False:
            volume_pair[a_pair[::-1]] = [volume, edvolume]
    
    else:
        volume_pair[a_pair]       = volume
        if polar == False:
            volume_pair[a_pair[::-1]] = volume
    
    return volume_pair


def pdb_dotarray_volume(contact_dict:dict, contacts_center:dict, contacts_dict:dict, polar:bool, weighted:bool, interval:float=0.2, disable_print=False):
    @timer(disable_print=disable_print)
    def count_volume():
        Pairs_Volume = []
        for a_pair in tqdm(contact_dict, disable=disable_print):
            # A pair of atoms, and all atoms in contact with them
            a_contact_dict = contact_dict[a_pair]
            a_pair_xyz  = {name: contacts_center[name] for name in a_contact_dict if name in contacts_center}   # Coordinates
            a_pair_info = {name: contacts_dict  [name] for name in a_contact_dict if name in contacts_dict  }   # Radius, volume
            
            volume_pair = determine_inner(a_pair, a_pair_xyz, a_pair_info, polar, weighted, interval)
            Pairs_Volume.append(volume_pair)
        
        volume = {} # Use update() to quickly merge dictionaries
        for pair in Pairs_Volume:
            volume.update(pair)

        return volume
    
    return count_volume()
