from Bio.PDB import *
from COM_helix import COM_helix
from COM import COM_protein
import os
import numpy as np
import argparse
import pandas as pd


def prot_hel_dist(pdb_file):
    hel_COMs = COM_helix(pdb_file)
    prot_COM = COM_protein(pdb_file)
    prot_hel_dists = []
    for i in range(0, len(hel_COMs)):
        vect = np.array(hel_COMs[i]) - np.array(prot_COM)
        prot_hel_dists.append(np.sqrt(np.sum(vect ** 2)))
    return prot_hel_dists


if __name__ == '__main__':
    dir = '/home/stephen/Desktop/PDB'  # Enter your PDB directory

    cols = ['prot_name']
    for elem in range(0, 12):
        cols.append('Helix_num_' + str(elem + 1))
    df = pd.DataFrame(columns=cols)

    for filename in os.listdir(dir):
        prothel = prot_hel_dist(os.path.join(dir, filename))

        data = [filename]
        for elem in prothel:
            data.append(elem)

        df = df.append(pd.Series(data, index=cols[0:len(data)]), ignore_index=True)
    print(df)