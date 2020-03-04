from Bio.PDB import *
from COM_helix import COM_helix
import numpy as np

def PairwiseSep(h1, h2, str_file):
    hel_COMs = COM_helix(str_file)
    vect = np.array(hel_COMs[h1-1]) - np.array(hel_COMs[h2-1])
    length = np.sqrt(np.sum(vect**2))
    return length
