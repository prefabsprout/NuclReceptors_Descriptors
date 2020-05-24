from Bio.PDB import *
from math import degrees
import pandas as pd

def ch_clamp_dist(pdb_file, charge_clamps):

    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    chain = structure[0]['A']

    clamp_vectors = {}

    for res in chain:
        if res.id[1] == charge_clamps[0] or res.id[1] == charge_clamps[1] or res.id[1] == charge_clamps[2]:
            for atom in res:
                if atom.get_name() == 'CA':
                    clamp_vectors[res.id[1]] = atom.get_vector()

    angles = {}
    for elem in range(len(clamp_vectors)):
        angles[f'{charge_clamps[elem]}-{charge_clamps[elem - 1]}-{charge_clamps[elem - 2]}'] = degrees(calc_angle(
            clamp_vectors[charge_clamps[elem]], clamp_vectors[charge_clamps[elem - 1]], clamp_vectors[charge_clamps[elem - 2]]))

    clamp_angles = pd.Series(angles).to_csv('clamp_angle_' + '1db1.pdb' + '.csv')

    return clamp_angles