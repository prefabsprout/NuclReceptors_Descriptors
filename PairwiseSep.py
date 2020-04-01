from Bio.PDB import *
from COM_helix import COM_helix
from os import path
import numpy as np
import argparse
import pandas as pd


def PairwiseSep(str_file):
    hel_COMs = COM_helix(str_file)
    pairwise_seps = []
    for i in range(1, len(hel_COMs)):
        vect = np.array(hel_COMs[i][0]) - np.array(hel_COMs[i - 1][0])
        pairwise_seps.append(np.sqrt(np.sum(vect ** 2)))
    return pairwise_seps


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='input_file',
                        required=True,
                        type=str)
    args = parser.parse_args()

    in_file_path = args.input_file

    PairwiseSep = PairwiseSep(in_file_path)
    prot_name = path.basename(in_file_path)

    cols = ['prot_name']
    for elem in enumerate(PairwiseSep):
        cols.append('Helix_num_' + str(elem[0] + 1) + '_PairwiseSep')

    data = [prot_name]
    for elem in PairwiseSep:
        data.append(elem)

    df = pd.DataFrame([data], columns=cols)
    df.to_csv('PairwiseSep.csv', mode='a', header=False)
