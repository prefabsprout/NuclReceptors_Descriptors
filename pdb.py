from Bio.PDB import *

# def acc_per_hel(pdb_file):

p = PDBParser()
structure = p.get_structure('1db1', '1db1.pdb')
model = structure[0]
dssp = DSSP(model, '1db1.pdb')

helix_borders = []

for i in range(0, len(dssp.keys())):
    if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
        border = [dssp[list(dssp.keys())[i]][0]]
    elif dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i + 1]][2] != 'H':
        border.append(dssp[list(dssp.keys())[i]][0])
        helix_borders.append(border)

helix_content = []

for i in range(0, len(helix_borders)):
    elements = [el for el in range(helix_borders[i][0], helix_borders[i][1] + 1)]
    helix_content.append(elements)


with open('1db1.dssp') as dssp_file:
    dssp_lines = [line.rstrip().split() for line in dssp_file]

    for i in range(28, len(dssp_lines)):
        for j in range(0, len(helix_content)):
            print(list(dssp_lines[i]))
            print(helix_content[j])
            if dssp_lines[i][0] in helix_content[j]:
                словарь_в_значение += dssp_lines[i][6]

    # тут надо как-то высчитать для каждой спирали
    # надо словарь сделать для каждой спирали



    # for i in range(0, len(helix_content)):
    #     for j in range(helix_content[i][0], helix_content[i][1]+1):
    #         print(j, dssp[list(dssp.keys())[j]][3])

# DSSP data is accessed by a tuple (chain_id, res_id)
# a_key = list(dssp.keys())[]
# print(a_key)
# (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
# NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
# NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)

#print(dssp[a_key])
