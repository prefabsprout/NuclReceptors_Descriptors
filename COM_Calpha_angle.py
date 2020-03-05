from Bio.PDB import *
from COM import COM_protein
import numpy as np


def COM_Calpha_angle(str_file, hel):
    protCOM = COM_protein(str_file)
    parser = MMCIFParser()
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

    OHel1 = np.array(protCOM) - CA_coords[0]
    OHel2 = np.array(protCOM) - CA_coords[1]

    OHels = np.prod(OHel1) + np.prod(OHel2)

    OHel1abs = np.sqrt(np.sum(OHel1 ** 2))
    OHel2abs = np.sqrt(np.sum(OHel2 ** 2))
    return OHels / (OHel1abs * OHel2abs)
