from Bio.PDB import *
import pandas as pd


def ch_clamp_dist(pdb_file, charge_clamps):
    """Calculation of distance between charge clamp residues"""
    
    # getting data from pdb-file
    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    chain = structure[0]['A']
    
    # extracting vectors of coordinates for every residue in charge clamp
    clamp_vectors = {}

    for res in chain:
        if res.id[1] == charge_clamps[0] or res.id[1] == charge_clamps[1] or res.id[1] == charge_clamps[2]:
            for atom in res:
                if atom.get_name() == 'CA':
                    clamp_vectors[res.id[1]] = atom.get_vector()
    table = {}
    for elem in range(len(clamp_vectors)):
        table[f'{charge_clamps[elem]}-{charge_clamps[elem-1]}'] = clamp_vectors[charge_clamps[elem]] - clamp_vectors[charge_clamps[elem-1]]
    
    # calculation of distance between charge clamp residues
    dist = {}

    for line in table:
        dist[line] = (table[line][0]**2 + table[line][1]**2 + table[line][2]**2) ** 0.5

    return dist
