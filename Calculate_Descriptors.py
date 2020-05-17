from itertools import chain
from math import degrees
from Bio.PDB import *
import numpy as np
import os
import pandas as pd


def COM_protein(pdb_file):
    ATOMIC_WEIGHTS = {'H': 1.008, 'HE': 4.002602, 'LI': 6.94, 'BE': 9.012182,
                      'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.9984032,
                      'NE': 20.1797, 'NA': 22.98976928, 'MG': 24.305, 'AL': 26.9815386,
                      'SI': 28.085, 'P': 30.973762, 'S': 32.06, 'CL': 35.45, 'AR': 39.948,
                      'K': 39.0983, 'CA': 40.078, 'SC': 44.955912, 'TI': 47.867, 'V': 50.9415,
                      'CR': 51.9961, 'MN': 54.938045, 'FE': 55.845, 'CO': 58.933195,
                      'NI': 58.6934, 'CU': 63.546, 'ZN': 65.38, 'GA': 69.723, 'GE': 72.630,
                      'AS': 74.92160, 'SE': 78.96, 'BR': 79.904, 'RB': 85.4678, 'SR': 87.62,
                      'Y': 88.90585, 'ZR': 91.224, 'NB': 92.90638, 'MO': 95.96, 'TC': 98,
                      'RU': 101.07, 'RH': 102.90550, 'PD': 106.42, 'AG': 107.8682, 'CD': 112.411,
                      'IN': 114.818, 'SN': 118.710, 'SB': 121.760, 'TE': 127.60, 'I': 126.90447,
                      'XE': 131.293, 'CS': 132.9054519, 'BA': 137.327, 'LA': 138.90547,
                      'CE': 140.116, 'PR': 140.90765, 'ND': 144.242, 'PM': 145, 'SM': 150.36,
                      'EU': 151.964, 'GD': 157.25, 'TB': 158.92535, 'DY': 162.500, 'HO': 164.93032,
                      'ER': 167.259, 'TM': 168.93421, 'YB': 173.054, 'LU': 174.9668, 'HF': 178.49,
                      'TA': 180.94788, 'W': 183.84, 'RE': 186.207, 'OS': 190.23, 'IR': 192.217,
                      'PT': 195.084, 'AU': 196.966569, 'HG': 200.592, 'TL': 204.38, 'PB': 207.2,
                      'BI': 208.98040, 'PO': 209, 'AT': 210, 'RN': 222, 'FR': 223, 'RA': 226,
                      'AC': 227, 'TH': 232.03806, 'PA': 231.03588, 'U': 238.02891, 'NP': 237,
                      'PU': 244, 'AM': 243, 'CM': 247, 'BK': 247, 'CF': 251, 'ES': 252, 'FM': 257,
                      'MD': 258, 'NO': 259, 'LR': 262, 'RF': 267, 'DB': 268, 'SG': 269, 'BH': 270,
                      'HS': 269, 'MT': 278, 'DS': 281, 'RG': 281, 'CN': 285, 'UUT': 286, 'FL': 289,
                      'UUP': 288, 'LV': 293, 'UUS': 294}

    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    atom_struct = structure.get_atoms()

    atoms = [([coord * ATOMIC_WEIGHTS[atom.get_name()[0]] for coord in list(atom.get_coord())]) for atom in atom_struct]

    atom_struct = structure.get_atoms()
    total_mass = sum([ATOMIC_WEIGHTS[atom.get_name()[0]] for atom in atom_struct])

    COM = [coord / total_mass for coord in np.sum(atoms, axis=0)]
    return COM


def COM_Calpha_angle(pdb_file):
    protCOM = COM_protein(pdb_file)
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    model = structure[0]
    dssp = DSSP(model, pdb_file)

    helix_content = []

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
            border = [dssp.keys()[i][1][1]]
        elif dssp[list(dssp.keys())[i]][2] != 'H' and dssp[list(dssp.keys())[i - 1]][2] == 'H':
            border.append(dssp.keys()[i - 1][1][1])
            helix_content.append(border)
    CA_coords = []

    for elem in helix_content:
        helices = [list(elem)]
        model = structure[0]
        chain = model['A']
        for elem in helices:
            coord = []
            for res in elem:
                residue = chain[res]
                for atom in residue.get_atoms():
                    if atom.get_name() == 'CA':
                        coord.append(atom.get_coord())
            CA_coords.append(coord)

    angles = []
    for elem in CA_coords:
        OHel1 = elem[0] - np.array(protCOM)
        OHel2 = elem[1] - np.array(protCOM)

        OHels = np.dot(OHel1, OHel2)

        OHel1abs = np.linalg.norm(OHel1)
        OHel2abs = np.linalg.norm(OHel2)
        angles.append(degrees(OHels / (OHel1abs * OHel2abs)))
    return angles


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


def COM_helix(pdb_file):
    ATOMIC_WEIGHTS = {'H': 1.008, 'HE': 4.002602, 'LI': 6.94, 'BE': 9.012182,
                      'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.9984032,
                      'NE': 20.1797, 'NA': 22.98976928, 'MG': 24.305, 'AL': 26.9815386,
                      'SI': 28.085, 'P': 30.973762, 'S': 32.06, 'CL': 35.45, 'AR': 39.948,
                      'K': 39.0983, 'CA': 40.078, 'SC': 44.955912, 'TI': 47.867, 'V': 50.9415,
                      'CR': 51.9961, 'MN': 54.938045, 'FE': 55.845, 'CO': 58.933195,
                      'NI': 58.6934, 'CU': 63.546, 'ZN': 65.38, 'GA': 69.723, 'GE': 72.630,
                      'AS': 74.92160, 'SE': 78.96, 'BR': 79.904, 'RB': 85.4678, 'SR': 87.62,
                      'Y': 88.90585, 'ZR': 91.224, 'NB': 92.90638, 'MO': 95.96, 'TC': 98,
                      'RU': 101.07, 'RH': 102.90550, 'PD': 106.42, 'AG': 107.8682, 'CD': 112.411,
                      'IN': 114.818, 'SN': 118.710, 'SB': 121.760, 'TE': 127.60, 'I': 126.90447,
                      'XE': 131.293, 'CS': 132.9054519, 'BA': 137.327, 'LA': 138.90547,
                      'CE': 140.116, 'PR': 140.90765, 'ND': 144.242, 'PM': 145, 'SM': 150.36,
                      'EU': 151.964, 'GD': 157.25, 'TB': 158.92535, 'DY': 162.500, 'HO': 164.93032,
                      'ER': 167.259, 'TM': 168.93421, 'YB': 173.054, 'LU': 174.9668, 'HF': 178.49,
                      'TA': 180.94788, 'W': 183.84, 'RE': 186.207, 'OS': 190.23, 'IR': 192.217,
                      'PT': 195.084, 'AU': 196.966569, 'HG': 200.592, 'TL': 204.38, 'PB': 207.2,
                      'BI': 208.98040, 'PO': 209, 'AT': 210, 'RN': 222, 'FR': 223, 'RA': 226,
                      'AC': 227, 'TH': 232.03806, 'PA': 231.03588, 'U': 238.02891, 'NP': 237,
                      'PU': 244, 'AM': 243, 'CM': 247, 'BK': 247, 'CF': 251, 'ES': 252, 'FM': 257,
                      'MD': 258, 'NO': 259, 'LR': 262, 'RF': 267, 'DB': 268, 'SG': 269, 'BH': 270,
                      'HS': 269, 'MT': 278, 'DS': 281, 'RG': 281, 'CN': 285, 'UUT': 286, 'FL': 289,
                      'UUP': 288, 'LV': 293, 'UUS': 294}

    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    model = structure[0]
    dssp = DSSP(model, pdb_file)

    helices = []

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
            border = [dssp.keys()[i][1][1]]
        elif dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] == 'H':
            border.append(dssp.keys()[i][1][1])
        elif dssp[list(dssp.keys())[i]][2] != 'H' and dssp[list(dssp.keys())[i - 1]][2] == 'H':
            helices.append(border)

    hel_COM = []

    for elem in helices:
        helix_content = [list(elem)]
        model = structure[0]
        chain = model['A']
        helix_mass = 0
        weighted_coord = list()
        helix_com = list()

        for elem in helix_content:
            for res in elem:
                residue = chain[res]

                for atom in residue.get_atoms():
                    weight = ATOMIC_WEIGHTS[atom.get_name()[0]]
                    helix_mass += weight
                    weighted_coord.append([coord * weight for coord in list(atom.get_coord())])

            helix_com.append([coord / helix_mass for coord in np.sum(weighted_coord, axis=0)])
            hel_COM.append(helix_com)

    return hel_COM


def PairwiseSep(pdb_file):
    hel_COMs = list(chain(*COM_helix(pdb_file)))
    pairwise_seps = []
    for j in range(0, len(hel_COMs) - 1):
        result = []
        for i in range(j + 1, len(hel_COMs)):
            vect = np.array(hel_COMs[i]) - np.array(hel_COMs[j])
            result.append(np.sqrt(np.sum(vect ** 2)))
        pairwise_seps.append(result)
    return pairwise_seps


def prot_hel_dist(pdb_file):
    hel_COMs = COM_helix(pdb_file)
    prot_COM = COM_protein(pdb_file)
    prot_hel_dists = []
    for i in range(0, len(hel_COMs)):
        vect = np.array(hel_COMs[i]) - np.array(prot_COM)
        prot_hel_dists.append(np.sqrt(np.sum(vect ** 2)))
    return prot_hel_dists


if __name__ == '__main__':
    dir = '/home/stepan/Desktop/pdb'  # Enter your PDB directory

    cols = ['prot_name', 'COM_x', 'COM_y', 'COM_z']
    for elem in range(0, 12):
        cols.append('Angle_between_COM_and_Calpha_of_hel_' + str(elem + 1))
    cols = ['prot_name', 'dist_between_COM_clamp_1', 'dist_between_COM_clamp_2', 'dist_between_COM_clamp_3']
    for elem in range(0, 12):
        cols.append('Helix_x_num_' + str(elem + 1))
        cols.append('Helix_y_num_' + str(elem + 1))
        cols.append('Helix_z_num_' + str(elem + 1))
    for i in range(0, 11):
        for j in range(i + 1, 12):
            cols.append('Pairwise_sep_between_hel_' + str(i + 1) + '_and_hel_' + str(j + 1))
    for elem in range(0, 12):
        cols.append('Helix_num_' + str(elem + 1))
    df = pd.DataFrame(columns=cols)

    for filename in os.listdir(dir):
        data = [filename]

        COM_prot = COM_protein(os.path.join(dir, filename))
        for coord in COM_prot:
            data.append(coord)

        alpha_angle = COM_Calpha_angle(os.path.join(dir, filename))
        for elem in alpha_angle:
            data.append(elem)

        clamps = COM_clamp(os.path.join(dir, filename), [256, 257, 258])
        for dist in clamps:
            data.append(dist)

        pairseps = PairwiseSep(os.path.join(dir, filename))
        for elem in pairseps:
            for dist in elem:
                data.append(dist)

        prothel = prot_hel_dist(os.path.join(dir, filename))
        for elem in prothel:
            data.append(elem)

        df = df.append(pd.Series(data, index=cols[0:len(data)]), ignore_index=True)
    df.to_csv("calc_results.csv")
