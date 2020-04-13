from Bio.PDB import *
from COM_helix import COM_helix
from itertools import chain
import os
import numpy as np
import argparse
import pandas as pd


def PairwiseSep(str_file):
    hel_COMs = list(chain(*COM_helix(str_file)))
    pairwise_seps = []
    for j in range(0, len(hel_COMs)-1):
        result = []
        for i in range(j+1, len(hel_COMs)):
            vect = np.array(hel_COMs[i]) - np.array(hel_COMs[j])
            result.append(np.sqrt(np.sum(vect ** 2)))
        pairwise_seps.append(result)
    return pairwise_seps


if __name__ == '__main__':
    dir = '/home/stephen/Desktop/PDB'  # Enter your PDB directory
    for filename in os.listdir(dir):
        print(PairwiseSep(os.path.join(dir, filename)))
