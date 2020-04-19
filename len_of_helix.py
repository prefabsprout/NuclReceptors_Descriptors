from Bio.PDB import *


def len_of_hel(pdb_file, dssp_file):

    with open(dssp_file) as dssp_file:
        dssp_lines = [line.rstrip().split() for line in dssp_file]

    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    # рассчитываем границы спиралей по  dssp
    helix_borders = []

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
            border = [dssp[list(dssp.keys())[i]][0]]
        elif dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i + 1]][2] != 'H':
            border.append(dssp[list(dssp.keys())[i]][0])
            helix_borders.append(border)

    # для каждой границы спирали ставим в соответствие три координаты СА
    res = dict()

    for i in range(len(helix_borders)):
        res[i] = [dssp_lines[helix_borders[i][0]+27][-3], dssp_lines[helix_borders[i][0]+27][-2], dssp_lines[helix_borders[i][0]+27][-1]], [dssp_lines[helix_borders[i][1]+27][-3], dssp_lines[helix_borders[i][1]+27][-2], dssp_lines[helix_borders[i][1]+27][-1]]

    # расчитываем длину получившегося отрезка по формуле ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)^0.5
    coord = dict()
    lens_of_helices = dict()

    for key, value in res.items():
        coord[key] = [(float(value[1][0]) - float(value[0][0]))**2]
        coord[key].append((float(value[1][1]) - float(value[0][1]))**2)
        coord[key].append((float(value[1][2]) - float(value[0][2]))**2)
        lens_of_helices[key] = (coord[key][0] + coord[key][1] + coord[key][2]) ** 0.5
    
    return lens_of_helices
