import warnings
warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser, PPBuilder

from .PDB_constants import AA
from .utils import timer
from .typing import *

def pdb_dihedral_angle(name:str, file:FileLike, disable_print=False) -> Angles:
    @timer(disable_print=disable_print)
    def count_dihedral_angle():
        # Parse PDB file
        parser = PDBParser()
        structure = parser.get_structure(name, file)

        # Get peptides from the structure
        Angles = {}
        for pp in PPBuilder().build_peptides(structure):
            phi_psi_list = pp.get_phi_psi_list()

            for i, residue in enumerate(pp):
                chain_id = residue.get_parent().get_id()    # Get chain id
                residue_id = residue.get_id()[1]            # Get residue number
                residue_name = residue.get_resname()        # Get residue name
                
                ID = chain_id + str(residue_id) + '-' + AA[residue_name]
                Angles[ID] = phi_psi_list[i]

        return Angles

    return count_dihedral_angle()

