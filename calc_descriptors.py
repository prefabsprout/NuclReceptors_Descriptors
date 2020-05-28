from math import degrees
from Bio.PDB import *
import os
import pandas as pd


def sseCalc(pdb_file):
    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    resamount = len(dssp.keys()) + 1
    dssp_structures = ['H', 'B', 'E', 'G', 'I', 'T', 'S']
    sses = list()

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] in dssp_structures:
            sses.append(dssp[list(dssp.keys())[i]][2])

    sse = {'Helix': sses.count('H') / resamount * 100,
           'Beta bridge': sses.count('B') / resamount * 100,
           'Strand': sses.count('E') / resamount * 100,
           'Helix-3': sses.count('G') / resamount * 100,
           'Helix-5': sses.count('I') / resamount * 100,
           'Turn': sses.count('T') / resamount * 100,
           'Bend': sses.count('S') / resamount * 100,
           'Other': (resamount - len(sses)) / resamount * 100}

    # sse_percent = pd.Series(sse).to_csv(pdb_file+'.csv')
    return sse


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


def len_of_helices(pdb_file):

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

    # lens_hel = pd.Series(lens_of_helices).to_csv('lengths_'+pdb_file+'.csv')

    return lens_of_helices


def cos_hel(pdb_file):

    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    helix_borders = []

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
            border = [dssp.keys()[i][1][1]]
        elif dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i + 1]][2] != 'H':
            border.append(dssp.keys()[i][1][1])
            helix_borders.append(border)

    helices = dict()
    for i in range(0, len(helix_borders)):
        helices[i] = [structure[0]['A'][res]['CA'].get_vector() for res in [helix_borders[i]][0]]

    for elem in helices:
        helices[elem] = helices[elem][0]-helices[elem][1]

    cos_between_hel = dict()
    for el in range(0, len(helices) - 1):
        for an_el in range(el + 1, len(helices)):
            cos_between_hel[f'{el}-{an_el}'] = (helices[el] * helices[an_el]) / \
                                               (((helices[el][0] ** 2 + helices[el][1] ** 2 + helices[el][
                                                   2] ** 2) ** 0.5) *
                                                ((helices[an_el][0] ** 2 + helices[an_el][1] ** 2 + helices[an_el][
                                                    2] ** 2) ** 0.5))
    # data = pd.DataFrame(cos_between_hel).to_csv('cos_'+pdb_file+'.csv')
    return cos_between_hel


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

    # clamp_dist = pd.Series(dist).to_csv('clamp_dist_'+pdb_file+'.csv')
    return dist



def ch_clamp_angles(pdb_file, charge_clamps):

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

    # clamp_angles = pd.Series(angles).to_csv('clamp_angle_' + '1db1.pdb' + '.csv')

    return angles


if __name__ == '__main__':
    def calc_all(directory, db_output, clamp_resid, species, prep):
        dir = directory  # Enter your PDB directory

        cols = ['prot_name', 'species_name', "preparation"]

        cols_sse = ['prot_name', 'Helix', 'Beta bridge', 'Strand', 'Helix-3', 'Helix-5', 'Turn', 'Bend', 'Other']
        cols_acc = ['prot_name']
        for elem in range(1, 13):
            cols_acc.append(f'Acc_per_hel_{elem}')
        cols_len = ['prot_name']
        for elem in range(1, 13):
            cols_len.append(f'len_of_hel_{elem}')
        cols_cos = ['prot_name']
        for i in range(1, 12):
            for j in range(i + 1, 13):
                cols_cos.append(f'cos_between_hel_{i}_and_hel_{j}')
        cols_cl_dist = ['prot_name']
        for elem in range(1, 4):
            cols_cl_dist.append(f'dist_clamp_{elem}')
        cols_cl_angle = ['prot_name']
        for elem in range(1,3):
            for el in range(elem+1, 4):
                cols_cl_angle.append(f'clamp_angle_{elem}-{el}')

        df = pd.DataFrame(columns=cols)
        df_cl_dist = pd.DataFrame(columns=cols_cl_dist)
        df_cl_angles = pd.DataFrame(columns=cols_cl_angle)
        df_sse = pd.DataFrame(columns=cols_sse)
        df_acc = pd.DataFrame(columns=cols_acc)
        df_len = pd.DataFrame(columns=cols_len)
        df_cos = pd.DataFrame(columns=cols_cos)

        for filename in os.listdir(dir):
            data = [filename, species, prep]
            df = df.append(pd.Series(data, index=cols[0:len(data)]), ignore_index=True)

            try:
                sse = sseCalc(os.path.join(dir, filename))
            except:
                KeyError
            data_sse = [filename]
            for struct in sse:
                data_sse.append(sse[struct])
            df_sse = df_sse.append(pd.Series(data_sse, index=cols_sse[0:len(data_sse)]), ignore_index=True)

            try:
                acc = acc_per_hel(os.path.join(dir, filename))
            except:
                KeyError
            data_acc = [filename]
            for elem in acc:
                data_acc.append(acc[elem])
            df_acc = df_acc.append(pd.Series(data_acc, index=cols_acc[0:len(data_acc)]), ignore_index=True)

            try:
                lens_hels = len_of_helices(os.path.join(dir, filename))
            except:
                KeyError
            data_lens = [filename]
            for lens in lens_hels:
                data_lens.append(lens_hels[lens])
            df_len = df_len.append(pd.Series(data_lens, index=cols_len[0:len(data_lens)]), ignore_index=True)

            try:
                cos = cos_hel(os.path.join(dir, filename))
            except:
                KeyError
            data_cos = [filename]
            for elem in cos:
                data_cos.append(cos[elem])
            df_cos = df_cos.append(pd.Series(data_cos, index=cols_cos[0:len(data_cos)]), ignore_index=True)

            try:
                clamp_dist = ch_clamp_dist(os.path.join(dir, filename), clamp_resid)
            except:
                KeyError
            cl_dist = [filename]
            if clamp_dist is not None:
                for elem in clamp_dist:
                    cl_dist.append(clamp_dist[elem])
            df_cl_dist = df_cl_dist.append(pd.Series(cl_dist, index=cols_cl_dist[0:len(cl_dist)]),
                                         ignore_index=True)

            try:
                clamp_angle = ch_clamp_angles(os.path.join(dir, filename), clamp_resid)
            except:
                KeyError
            cl_angle = [filename]
            if clamp_angle is not None:
                for elem in clamp_angle:
                    cl_angle.append(clamp_angle[elem])
            df_cl_angles = df_cl_angles.append(pd.Series(cl_angle, index=cols_cl_angle[0:len(cl_angle)]),
                                           ignore_index=True)

            df_concat = df.merge(df_sse, on='prot_name').merge(df_acc, on='prot_name') \
                .merge(df_len, on='prot_name').merge(df_cos, on='prot_name') \
                .merge(df_cl_dist, on='prot_name').merge(df_cl_angles, on='prot_name')
        df_concat.to_csv(db_output)

    calc_all("/home/masha/Descriptors/VDR_PDB/Danio rerio", "calc_results_danio.csv", [274, 292, 446],
             "zebrafish", 0)
    calc_all("/home/masha/Descriptors/VDR_PDB/Homo sapiens", "calc_results_homosap.csv", [246, 264, 420],
             "human", 0)
    calc_all("/home/masha/Descriptors/VDR_PDB/Rattus norvegicus", "calc_results_rattus.csv", [242, 260, 416],
             "rat", 0)

    calc_all("/home/masha/Descriptors/VDR_PDB_prep/Danio rerio", "calc_results_danio_prep.csv", [274, 292, 446],
             "zebrafish", 1)
    calc_all("/home/masha/Descriptors/VDR_PDB_prep/Homo sapiens", "calc_results_homosap_prep.csv",
             [246, 264, 420],
             "human", 1)
    calc_all("/home/masha/Descriptors/VDR_PDB_prep/Rattus norvegicus", "calc_results_rattus_prep.csv",
             [242, 260, 416], "rat", 1)
