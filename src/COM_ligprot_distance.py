from COM import COM_protein
from COM_ligand import *
import numpy as np

def prot_ligand_dist(pdb_file,lig_file):
    """Caclulate distance between protein's center of mass and ligand's center of mass"""
    # Calculate protein's center of mass and ligand's center of mass
    hel_COMs = COM_protein(pdb_file)
    ligand_COM = COM_ligand(pdb_file)
    # Calculte protein's COM and ligand's center of mass as vector
    vect = np.array(hel_COMs) - np.array(ligand_COM)
    return np.sqrt(np.sum(vect ** 2))
