from Bio.PDB import *
from COM import COM_protein
import os
import numpy as np
import pandas as pd
import argparse


def COM_clamp(pdb_file, ch_clamps):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    model = structure[0]
    chain = model['A']

    atom_coord = []
    com = np.array(COM_protein(pdb_file))

    for res in chain:
        if res.id[1] == ch_clamps[0] or res.id[1] == ch_clamps[1] or res.id[1] == ch_clamps[2]:
            for atom in res:
                if atom.get_name() == 'CA':
                    atom_coord.append(atom.get_coord())

    ch_clamp_dist = []
    for elem in atom_coord:
        vect = np.array(elem) - np.array(com)
        ch_clamp_dist.append(np.sqrt(np.sum(vect ** 2)))

    return ch_clamp_dist


if __name__ == '__main__':
    dir = '/home/stephen/Desktop/PDB'  # Enter your PDB directory

    cols = ['prot_name', 'dist_between_COM_clamp_1', 'dist_between_COM_clamp_2', 'dist_between_COM_clamp_3']
    df = pd.DataFrame(columns=cols)

    for filename in os.listdir(dir):
        clamps = COM_clamp(os.path.join(dir, filename), [256,257,258])

        data = [filename]
        for dist in clamps:
            data.append(dist)

        df = df.append(pd.Series(data, index=cols[0:len(data)]), ignore_index=True)
        print(df)
