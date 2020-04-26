from Bio.PDB import *
from COM import COM_protein
from math import degrees
import os
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
    dir = '/home/stephen/Desktop/PDB'  # Enter your PDB directory

    cols = ['prot_name']
    for elem in range(0, 12):
        cols.append('Angle_between_COM_and_Calpha_of_hel_' + str(elem + 1))
    df = pd.DataFrame(columns=cols)

    for filename in os.listdir(dir):
        alpha_angle = COM_Calpha_angle(os.path.join(dir, filename))

        data = [filename]
        for elem in alpha_angle:
            data.append(elem)

        df = df.append(pd.Series(data, index=cols[0:len(data)]), ignore_index=True)
