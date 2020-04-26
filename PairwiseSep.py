from Bio.PDB import *
from COM_helix import COM_helix
from itertools import chain
import os
import numpy as np
import argparse
import pandas as pd


def PairwiseSep(pdb_file):
    hel_COMs = list(chain(*COM_helix(pdb_file)))
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

    cols = ['prot_name']
    for i in range(0, 11):
        for j in range(i+1, 12):
            cols.append('Pairwise_sep_between_hel_' + str(i + 1) + '_and_hel_' + str(j+1))
    print(cols)

    df = pd.DataFrame(columns=cols)

    for filename in os.listdir(dir):
        pairseps = PairwiseSep(os.path.join(dir, filename))

        data = [filename]
        for elem in pairseps:
            for dist in elem:
                data.append(dist)

        df = df.append(pd.Series(data, index=cols[0:len(data)]), ignore_index=True)
