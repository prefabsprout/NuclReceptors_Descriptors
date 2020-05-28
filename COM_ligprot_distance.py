from COM import COM_protein
from COM_ligand import *
import numpy as np

def prot_hel_dist(pdb_file,lig_file):
    hel_COMs = COM_protein(pdb_file)
    prot_COM = COM_ligand(pdb_file)
    vect = np.array(hel_COMs[i]) - np.array(prot_COM)
    return np.sqrt(np.sum(vect ** 2))
