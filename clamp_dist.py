from Bio.PDB import *
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
    table = {}
    for elem in range(len(clamp_vectors)):
        table[f'{charge_clamps[elem]}-{charge_clamps[elem-1]}'] = clamp_vectors[charge_clamps[elem]] - clamp_vectors[charge_clamps[elem-1]]

    dist = {}

    for line in table:
        dist[line] = (table[line][0]**2 + table[line][1]**2 + table[line][2]**2) ** 0.5

    clamp_dist = pd.Series(dist).to_csv('clamp_dist_'+pdb_file+'.csv')

    return clamp_dist 