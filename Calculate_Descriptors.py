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
    helices = []
    for i in range(0, len(dssp.keys())):
        if dssp.keys()[i][0] == 'A':
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
    try:
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
    except:
        KeyError


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
        if dssp.keys()[i][0] == 'A':
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
    def calc_all(directory, db_output, clamp_resid, species, prep):
        dir = directory # Enter your PDB directory

        cols = ['prot_name', 'species_name', "preparation"]

        com_cols = ['prot_name','COM_x', 'COM_y', 'COM_z']

        angle_cols = ['prot_name']
        for elem in range(0, 12):
            angle_cols.append('Angle_between_COM_and_Calpha_of_hel_' + str(elem + 1))

        comclampdist_cols = ['prot_name','dist_between_COM_clamp_1', 'dist_between_COM_clamp_2', 'dist_between_COM_clamp_3']

        comhel_cols = ['prot_name']
        for elem in range(0, 12):
            comhel_cols.append('COM_helix_x_num_' + str(elem + 1))
            comhel_cols.append('COM_helix_y_num_' + str(elem + 1))
            comhel_cols.append('COM_helix_z_num_' + str(elem + 1))

        pairwise_cols = ['prot_name']
        for i in range(0, 11):
            for j in range(i + 1, 12):
                pairwise_cols.append('Pairwise_sep_between_hel_' + str(i + 1) + '_and_hel_' + str(j + 1))

        protheldist_cols = ['prot_name']
        for elem in range(0, 12):
            protheldist_cols.append('Prot_Helix_' + str(elem + 1) + '_distance')

        df = pd.DataFrame(columns=cols)
        df_clamps = pd.DataFrame(columns=comclampdist_cols)
        df_com = pd.DataFrame(columns=com_cols)
        df_alphaangle = pd.DataFrame(columns=angle_cols)
        df_comhels = pd.DataFrame(columns=comhel_cols)
        df_pairseps = pd.DataFrame(columns=pairwise_cols)
        df_prothel = pd.DataFrame(columns=protheldist_cols)

        for filename in os.listdir(dir):
            data = [filename, species, prep]
            df = df.append(pd.Series(data, index=cols[0:len(data)]), ignore_index=True)

            try:
                COM_prot = COM_protein(os.path.join(dir, filename))
            except:
                KeyError
            data_com = [filename]
            for coord in COM_prot:
                data_com.append(coord)
            df_com = df_com.append(pd.Series(data_com, index=com_cols[0:len(data_com)]), ignore_index=True)

            try:
                alpha_angle = COM_Calpha_angle(os.path.join(dir, filename))
            except:
                KeyError
            data_alphaagnle = [filename]
            for elem in alpha_angle:
                data_alphaagnle.append(elem)
            df_alphaangle = df_alphaangle.append(pd.Series(data_alphaagnle, index=angle_cols[0:len(data_alphaagnle)]), ignore_index=True)

            try:
                clamps = COM_clamp(os.path.join(dir, filename), clamp_resid)
            except:
                KeyError
            data_clamps = [filename]
            if clamps is not None:
                for dist in clamps:
                    data_clamps.append(dist)
            df_clamps = df_clamps.append(pd.Series(data_clamps, index=comclampdist_cols[0:len(data_clamps)]), ignore_index=True)

            try:
                com_hels = COM_helix(os.path.join(dir, filename))
            except:
                KeyError
            data_comhels = [filename]
            for elem in com_hels:
                for coord in elem[0]:
                    data_comhels.append(coord)
            df_comhels = df_comhels.append(pd.Series(data_comhels, index=comhel_cols[0:len(data_comhels)]), ignore_index=True)

            try:
                pairseps = PairwiseSep(os.path.join(dir, filename))
            except:
                KeyError
            data_pairseps = [filename]
            for elem in pairseps:
                for dist in elem:
                    data_pairseps.append(dist)
            df_pairseps = df_pairseps.append(pd.Series(data_pairseps, index=pairwise_cols[0:len(data_pairseps)]), ignore_index=True)

            try:
                prothel = prot_hel_dist(os.path.join(dir, filename))
            except:
                KeyError
            data_prothel = [filename]
            for elem in prothel:
                data_prothel.append(elem)
            df_prothel = df_prothel.append(pd.Series(data_prothel, index=protheldist_cols[0:len(data_prothel)]), ignore_index=True)

            df_concat = df.merge(df_prothel, on='prot_name').merge(df_pairseps,on='prot_name')\
                .merge(df_comhels,on='prot_name').merge(df_clamps,on='prot_name')\
                .merge(df_alphaangle,on='prot_name').merge(df_com,on='prot_name')
        df_concat.to_csv(db_output)


    calc_all("/home/stepan/Git_Repo/VDR_PDB/Danio rerio", "calc_results_danio.csv", [274, 292, 446],
             "zebrafish",0)
    calc_all("/home/stepan/Git_Repo/VDR_PDB/Homo sapiens", "calc_results_homosap.csv", [246, 264, 420],
             "human",0)
    calc_all("/home/stepan/Git_Repo/VDR_PDB/Rattus norvegicus", "calc_results_rattus.csv", [242, 260, 416],
             "rat",0)

    calc_all("/home/stepan/Git_Repo/VDR_PDB_prep/Danio rerio", "calc_results_danio_prep.csv",[274, 292, 446],
             "zebrafish",1)
    calc_all("/home/stepan/Git_Repo/VDR_PDB_prep/Homo sapiens", "calc_results_homosap_prep.csv",[246, 264, 420],
             "human",1)
    calc_all("/home/stepan/Git_Repo/VDR_PDB_prep/Rattus norvegicus", "calc_results_rattus_prep.csv",[242, 260, 416],
             "rat",1)