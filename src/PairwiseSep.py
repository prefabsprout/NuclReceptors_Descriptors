from COM_helix import COM_helix
from itertools import chain
import numpy as np

def PairwiseSep(pdb_file):
    """Calculate separation distance between every helix"""
    # Calculate centers of mass of every helix
    hel_COMs = list(chain(*COM_helix(pdb_file)))
    # Calculate pairwise separations between every COM's of every helix as vector distance
    pairwise_seps = []
    for j in range(0, len(hel_COMs) - 1):
        result = []
        for i in range(j + 1, len(hel_COMs)):
            vect = np.array(hel_COMs[i]) - np.array(hel_COMs[j])
            result.append(np.sqrt(np.sum(vect ** 2)))
        pairwise_seps.append(result)
    return pairwise_seps