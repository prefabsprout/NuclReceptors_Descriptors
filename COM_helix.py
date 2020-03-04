from Bio.PDB import *
import numpy as np

def COM_helix(str_file):

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

    parser = MMCIFParser()
    structure = parser.get_structure('protein', str_file)

    model = structure[0]
    dssp = DSSP(model, str_file)

    helix_content = []
    helix_com = []

    for i in range(0,len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i-1]][2] != 'H':
            border = [dssp[list(dssp.keys())[i]][0]]
        elif dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i-1]][2] == 'H':
            border.append(dssp[list(dssp.keys())[i]][0])
        elif dssp[list(dssp.keys())[i]][2] != 'H' and dssp[list(dssp.keys())[i - 1]][2] == 'H':
            helix_content.append(border)

    for elem in helix_content:
        for residue in enumerate(structure.get_residues()):
            if (residue[0]+1) in elem:
               atom_struct = residue[1].get_atoms()
               atoms = [([coord * ATOMIC_WEIGHTS[atom.get_name()[0]] for coord in list(atom.get_coord())]) for atom in atom_struct]

               atom_struct = residue[1].get_atoms()
               total_mass = sum([ATOMIC_WEIGHTS[atom.get_name()[0]] for atom in atom_struct])
        helix_com.append([coord / total_mass for coord in np.sum(atoms, axis=0)])

    return helix_com


print(COM_helix('/home/stephen/Desktop/3b0t.cif'))