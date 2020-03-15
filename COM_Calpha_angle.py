from Bio.PDB import *
from COM import COM_protein
import numpy as np
import os


def COM_Calpha_angle(str_file, hel):
    protCOM = COM_protein(str_file)
    parser = PDBParser()
    structure = parser.get_structure('protein', str_file)

    model = structure[0]
    dssp = DSSP(model, str_file)

    helix_content = []

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
            border = [dssp[list(dssp.keys())[i]][0]]
        elif dssp[list(dssp.keys())[i]][2] != 'H' and dssp[list(dssp.keys())[i - 1]][2] == 'H':
            border.append(dssp[list(dssp.keys())[i]][0])
            helix_content.append(border)

    CA_coords = []

    for residue in enumerate(structure.get_residues()):
        if (residue[0] + 1) in helix_content[hel - 1]:
            atom_struct = residue[1].get_atoms()
            for atom in atom_struct:
                if atom.get_name() == 'CA':
                    CA_coords.append(atom.get_coord())

    OHel1 = CA_coords[0] - np.array(protCOM)
    OHel2 = CA_coords[1] - np.array(protCOM)

    OHels = np.dot(OHel1,OHel2)

    OHel1abs = np.linalg.norm(OHel1)
    OHel2abs = np.linalg.norm(OHel2)
    return OHels / (OHel1abs * OHel2abs)


dir = '/home/stephen/Desktop/PDB' # Enter your PDB directory
helnum = int() # Enter helix number
for filename in os.listdir(dir):
    print(COM_Calpha_angle(os.path.join(dir, filename),helnum))