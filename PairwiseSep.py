from Bio.PDB import *
from COM_helix import COM_helix
import numpy as np

def PairwiseSep(h1, h2, str_file):
    hel_COMs = COM_helix(str_file)
    vect = np.array(hel_COMs[h1-1]) - np.array(hel_COMs[h2-1])
    length = np.sqrt(np.sum(vect**2))
    return length

dir = '/home/stephen/Desktop/PDB' # Enter your PDB directory
h1 = int()
h2 = int()
for filename in os.listdir(dir):
    print(PairwiseSep(h1,h2, os.path.join(dir, filename)))