from Bio.PDB import *
from COM import COM_protein
from math import degrees
import numpy as np

def COM_Calpha_angle(pdb_file):
    """Calculate angles between protein's center of mass and alpha carbon atom of every helix"""
    # Initialize PDB structure
    protCOM = COM_protein(pdb_file)
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    # Initailize DSSP algorithm and append it to PDB structure
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    # Form list of helix (every nested list represents a helix and amino acids it contains)
    helix_content = []
    for i in range(0, len(dssp.keys())):
        # We need alpha helices. In DSSP output there is a column which represents a type of
        # secondary structure. So we need helix - H.
        if dssp.keys()[i][0] == 'A':
            if dssp[list(dssp.keys())[i]][2] == 'H' and dssp[list(dssp.keys())[i - 1]][2] != 'H':
                border = [dssp.keys()[i][1][1]]
            elif dssp[list(dssp.keys())[i]][2] != 'H' and dssp[list(dssp.keys())[i - 1]][2] == 'H':
                border.append(dssp.keys()[i-1][1][1])
                helix_content.append(border)

    # Calculate coordinates of alpha carbon atom of every helix
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

    # Calculate angles between protein's center of mass and helices of every alpha carbon
    angles = []
    for elem in CA_coords:
        # Calculate it as angle formed by two vectors
        OHel1 = elem[0] - np.array(protCOM)
        OHel2 = elem[1] - np.array(protCOM)

        OHels = np.dot(OHel1,OHel2)

        OHel1abs = np.linalg.norm(OHel1)
        OHel2abs = np.linalg.norm(OHel2)
        angles.append(degrees(OHels / (OHel1abs * OHel2abs)))
    return angles
