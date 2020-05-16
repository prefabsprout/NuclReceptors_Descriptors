from Bio.PDB import *
import pandas as pd


def len_of_hel(pdb_file):

    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    # рассчитываем границы спиралей по  dssp
    helix_borders = []

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
            border = [dssp.keys()[i][1][1]]
        elif dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i + 1]][2] != 'H':
            border.append(dssp.keys()[i][1][1])
            helix_borders.append(border)

    # для каждой граничной точки экстрагируем вектор координат
    helices = {}
    for i in range(0, len(helix_borders)):
        helices[i] = [structure[0]['A'][res]['CA'].get_vector() for res in [helix_borders[i]][0]]

    for elem in helices:
        helices[elem] = (helices[elem][1]-helices[elem][0])

    lens_of_helices = dict()

    # расчитываем длину получившегося отрезка по формуле ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)^0.5
    for el in helices:
        lens_of_helices[el] = (helices[el][0]**2 + helices[el][1]**2 + helices[el][2]**2) ** 0.5

    lens_hel = pd.Series(lens_of_helices).to_csv(pdb_file+'.csv')

    return lens_hel 
