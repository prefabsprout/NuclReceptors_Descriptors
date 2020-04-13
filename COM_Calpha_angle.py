from Bio.PDB import *
from COM import COM_protein
from math import degrees
from os import path
import numpy as np
import pandas as pd
import argparse

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
            border.append(dssp.keys()[i-1][1][1])
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

        OHels = np.dot(OHel1,OHel2)

        OHel1abs = np.linalg.norm(OHel1)
        OHel2abs = np.linalg.norm(OHel2)
        angles.append(degrees(OHels / (OHel1abs * OHel2abs)))
    return angles


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='input_file',
                        required=True,
                        type=str)
    args = parser.parse_args()

    in_file_path = args.input_file

    COM_Сalpha = COM_Calpha_angle(in_file_path)
    prot_name = path.basename(in_file_path)

    cols = ['prot_name']
    for elem in enumerate(COM_Сalpha):
        cols.append('Helix_num_' + str(elem[0] + 1) + '_COM_Calpha_angle')

    data = [prot_name]
    for elem in COM_Сalpha:
        data.append(elem)

    df = pd.DataFrame([data], columns=cols)
    df.to_csv('COM_Calpha_angle.csv', mode='a', header=False)
