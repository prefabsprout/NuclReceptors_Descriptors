from Bio.PDB import *
from COM_helix import COM_helix
import numpy as np

hel_COMs = COM_helix('/home/stephen/Desktop/3b0t.cif')

def PairwiseSep(h1, h2):
    vect = np.sqrt((np.array(hel_COMs[h1-1]) - np.array(hel_COMs[h2-1]))**2)
    length = np.sqrt(np.sum(vect**2))
    return length
