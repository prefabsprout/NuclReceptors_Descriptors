from Bio.PDB import *
from COM import COM_protein
import numpy as np


def COM_clamp(pdb_file, ch_clamps):
    """Calculate distances between protein's center of mass and every charge clamps"""
    try:
        # Initialize PDB structure
        parser = PDBParser()
        structure = parser.get_structure('protein', pdb_file)

        model = structure[0]
        chain = model['A']

        atom_coord = []
        com = np.array(COM_protein(pdb_file)) # We need here to calculate protein's center of mass

        # Find charge clamps from user input in PDB structure and get it coordinates
        for res in chain:
            if res.id[1] == ch_clamps[0] or res.id[1] == ch_clamps[1] or res.id[1] == ch_clamps[2]:
                for atom in res:
                    if atom.get_name() == 'CA':
                        atom_coord.append(atom.get_coord())

        # Calculate distance between center of mass and every charge clamps
        ch_clamp_dist = []
        for elem in atom_coord:
            vect = np.array(elem) - np.array(com)
            ch_clamp_dist.append(np.sqrt(np.sum(vect ** 2)))

        return ch_clamp_dist
    except:
        KeyError