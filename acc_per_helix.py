from Bio.PDB import *

def acc_per_hel(pdb_file):

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
        helices[i] = {el : dssp[list(dssp.keys())[el-1]][3] for el in range(helix_borders[i][0], helix_borders[i][1]+1)}

    for key in helices.keys():
        for res in helices[key].keys():
            if key in acc:
                acc[key] += helices[key][res] * (len(dssp.keys())+1)
            else:
                acc[key] = helices[key][res] * (len(dssp.keys())+1)

    for key, acc_sum in acc.items():
        acc[key] = acc_sum/(len(helices[key]))

    return acc