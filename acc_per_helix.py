from Bio.PDB import *


def acc_per_hel(pdb_file):
    res_max_acc = {
        "A": 106.0,
        "R": 248.0,
        "N": 157.0,
        "D": 163.0,
        "C": 135.0,
        "Q": 198.0,
        "E": 194.0,
        "G": 84.0,
        "H": 184.0,
        "I": 169.0,
        "L": 164.0,
        "K": 205.0,
        "M": 188.0,
        "F": 197.0,
        "P": 136.0,
        "S": 130.0,
        "T": 142.0,
        "W": 227.0,
        "Y": 222.0,
        "V": 142.0,
    }

    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    helix_borders = []

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
            border = [dssp[list(dssp.keys())[i]][0]]
        elif dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i + 1]][2] != 'H':
            border.append(dssp[list(dssp.keys())[i]][0])
            helix_borders.append(border)

    helices = dict()
    acc = dict()

    for i in range(len(helix_borders)):
        helices[i] = {el : dssp[list(dssp.keys())[el-1]][3] * res_max_acc[dssp[list(dssp.keys())[el-1]][1]] for el in range(helix_borders[i][0], helix_borders[i][1]+1)}

    for key in helices.keys():
        for res in helices[key].keys():
            if key in acc:
                acc[key] += helices[key][res]
            else:
                acc[key] = helices[key][res]

    for key, acc_sum in acc.items():
        acc[key] = acc_sum/(len(helices[key]))
    
    return acc
