from Bio.PDB import *
import pandas as pd


def sseCalc(pdb_file):
    """Calculation of secondary structure content"""
    p = PDBParser()
    structure = p.get_structure('protein', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)
    
    # preparing for extruction sse from dssp structure
    resamount = len(dssp.keys()) + 1
    dssp_structures = ['H', 'B', 'E', 'G', 'I', 'T', 'S']
    sses = list()
    
    # extracting sse from dssp
    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] in dssp_structures:
            sses.append(dssp[list(dssp.keys())[i]][2])
    
    # making dict of all possible secondary structures and counting their percentage
    sse = {'Helix': sses.count('H') / resamount * 100,
           'Beta bridge': sses.count('B') / resamount * 100,
           'Strand': sses.count('E') / resamount * 100,
           'Helix-3': sses.count('G') / resamount * 100,
           'Helix-5': sses.count('I') / resamount * 100,
           'Turn': sses.count('T') / resamount * 100,
           'Bend': sses.count('S') / resamount * 100,
           'Other': (resamount - len(sses)) / resamount * 100}

    return sse
