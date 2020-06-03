from Bio.PDB import *
from COM_helix import COM_helix
from COM import COM_protein
import numpy as np


def prot_hel_dist(pdb_file):
    """Calculate distance between protein's center of mass and between every helix's center of mass"""
    # Calculate protein's center of mass and centers of helices masses
    hel_COMs = COM_helix(pdb_file)
    prot_COM = COM_protein(pdb_file)
    # Calculate distance between protein's COM and helices COMs as vector distance
    prot_hel_dists = []
    for i in range(0, len(hel_COMs)):
        vect = np.array(hel_COMs[i]) - np.array(prot_COM)
        prot_hel_dists.append(np.sqrt(np.sum(vect ** 2)))
    return prot_hel_dists
