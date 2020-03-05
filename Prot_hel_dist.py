from Bio.PDB import *
from COM_helix import COM_helix
from COM import COM_protein
import numpy as np


def prot_hel_dist(h, str_file):
    hel_COMs = COM_helix(str_file)
    prot_COM = COM_protein(str_file)
    vect = np.array(hel_COMs[h - 1]) - np.array(prot_COM)
    length = np.sqrt(np.sum(vect ** 2))
    return length
