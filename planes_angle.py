from Bio.PDB import *
from sympy import *
import numpy as np
import argparse


def plane_angle(str_file, l1=[0, 1, 2], l2=[3, 4, 5], l3=[7, 8, 9]):
    parser = PDBParser()
    structure = parser.get_structure('protein', str_file)

    model = structure[0]
    dssp = DSSP(model, str_file)

    helices = []

    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] != 'H' and dssp[list(dssp.keys())[i - 1]][2] == 'H':
            helices.append([dssp.keys()[i][1][1]])

    last_atoms = []
    for elem in helices:
        helix_content = [list(elem)]
        model = structure[0]
        chain = model['A']
        helix_mass = 0
        coord = list()
        helix_com = list()

        for elem in helix_content:
            for res in elem:
                residue = chain[res]
                for atom in residue.get_atoms():
                    coord.append([coord for coord in list(atom.get_coord())])
                last_atoms.append(coord[-1])

    args = []
    for i in range(l1[0], l1[-1]+1):
        args.append(Point3D(last_atoms[i]))
    first_layer = Plane(*args)

    args = []
    for i in range(l2[0], l2[-1]+1):
        args.append(Point3D(last_atoms[i]))
    second_layer = Plane(*args)

    args = []
    for i in range(l3[0], l3[-1]+1):
        args.append(Point3D(last_atoms[i]))
    third_layer = Plane(*args)

    print('Angle between 1st and 2d layer:', N(first_layer.angle_between(second_layer)))
    print('Angle between 1st and 3d layer:', N(first_layer.angle_between(third_layer)))
    print('Angle between 2d and 3d layer:', N(second_layer.angle_between(third_layer)))

parser = argparse.ArgumentParser()

parser.add_argument('-i', dest='input_file',
                    required=True,
                    type=str)
args = parser.parse_args()

in_file_path = args.input_file
print(plane_angle(in_file_path))
