from Bio.PDB import *
from COM_helix import COM_helix
from COM import COM_protein
from os import path
import numpy as np
import argparse
import pandas as pd


def prot_hel_dist(str_file):
    hel_COMs = COM_helix(str_file)
    prot_COM = COM_protein(str_file)
    prot_hel_dists = []
    for i in range(1, len(hel_COMs)):
        vect = np.array(hel_COMs[i]) - np.array(prot_COM)
        prot_hel_dists.append(np.sqrt(np.sum(vect ** 2)))
    return prot_hel_dists


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='input_file',
                        required=True,
                        type=str)
    args = parser.parse_args()

    in_file_path = args.input_file

    prothel = prot_hel_dist(in_file_path)
    prot_name = path.basename(in_file_path)

    cols = ['prot_name']
    for elem in enumerate(prothel):
        cols.append('Helix_num_' + str(elem[0] + 1) + '_ProtHelDist')

    data = [prot_name]
    for elem in prothel:
        data.append(elem)

    df = pd.DataFrame([data], columns=cols)
    df.to_csv('Prot_hel_dist.csv', mode='a', header=False)
