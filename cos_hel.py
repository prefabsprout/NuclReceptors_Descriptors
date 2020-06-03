from Bio.PDB import *
import pandas as pd


def cos_hel(pdb_file):
    """Calculation of cos between all helices in structure"""
    
    # getting structure from pdb-file
    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)
    
    # extracting borders of helices from dssp and vectors for them
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
    
    # making vectors for helices
    for elem in helices:
        helices[elem] = helices[elem][0]-helices[elem][1]
    
    # calculation cos for every pair of helices in structure
    cos_between_hel = dict()
    for el in helices:
        cos_between_hel[el] = [(helices[el]*helices[an_el]) / \
                                (((helices[el][0]**2+helices[el][1]**2+helices[el][2]**2)**0.5) *
                                ((helices[an_el][0]**2+helices[an_el][1]**2+helices[an_el][2]**2)**0.5))
                                for an_el in helices]
    data = pd.DataFrame(cos_between_hel).to_csv('cos_'+pdb_file+'.csv')
    return data
