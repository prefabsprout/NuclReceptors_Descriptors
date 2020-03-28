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
            for atom in enumerate(res):
                if atom[0] == 0:
                    atom_coord.append(atom[1].get_coord())

    ch_clamp_dist = []
    for elem in atom_coord:
        vect = np.array(elem) - np.array(com)
        length = np.sqrt(np.sum(vect ** 2))
        ch_clamp_dist.append(length)
        
    return ch_clamp_dist

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', dest='input_file',
                        required=True,
                        type=str)
    
    parser.add_argument('--chargeclamps', dest='chcl',
                        nargs='+',
                        required=False,
                        type=int)
    
    args = parser.parse_args()
    
    in_file_path = args.input_file
    ch_clamps = args.chcl

    dists = COM_clamp(in_file_path, ch_clamps)
    prot_name = path.basename(in_file_path)

    data = {'protein_name': [prot_name],
            'clamp1_id': [ch_clamps[0]],
            'clamp2_id': [ch_clamps[1]],
            'clamp3_id': [ch_clamps[2]],
            'clamp1_dist': [dists[0]],
            'clamp1_dist': [dists[1]],
            'clamp1_dist': [dists[2]]}
    my_data = pd.DataFrame(data)
    my_data.to_csv('my_csv.csv', mode='a', header=False)
