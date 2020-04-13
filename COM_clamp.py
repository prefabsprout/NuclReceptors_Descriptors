from Bio.PDB import *
from COM import COM_protein
from os import path
import numpy as np
import pandas as pd
import argparse


def COM_clamp(str_file, ch_clamps):
    parser = PDBParser()
    structure = parser.get_structure('protein', str_file)

    model = structure[0]
    chain = model['A']

    atom_coord = []
    com = np.array(COM_protein(str_file))

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

print(COM_clamp('/home/stephen/Desktop/PDB/1db1.pdb', [256,257,258]))

#
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument('-i', dest='input_file',
#                         required=True,
#                         type=str)
#
#     parser.add_argument('--chargeclamps', dest='chcl',
#                         nargs='+',
#                         required=False,
#                         type=int)
#
#     args = parser.parse_args()
#
#     in_file_path = args.input_file
#     ch_clamps = args.chcl
#
#     dists = COM_clamp(in_file_path, ch_clamps)
#     prot_name = path.basename(in_file_path)
#
#     cols = ['prot_name']
#     for elem in enumerate(dists):
#         cols.append('Helix_num_' + str(elem[0] + 1) + '_COM_clamp')
#
#     data = [prot_name, ch_clamps[1], ch_clamps[2], ch_clamps[3]]
#     for elem in dists:
#         data.append(elem[0])
#         data.append(elem[1])
#         data.append(elem[2])
#
#     df = pd.DataFrame([data], columns=cols)
#     df.to_csv('COM_clamp.csv', mode='a', header=False)
