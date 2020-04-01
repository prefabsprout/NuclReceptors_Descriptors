from Bio.PDB import *
from sympy import *
from os import path
from math import degrees
import numpy as np
import argparse
import pandas as pd


def plane_angle(str_file, l1, l2, l3):
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
        coord = list()

        for elem in helix_content:
            for res in elem:
                residue = chain[res]
                for atom in residue.get_atoms():
                    coord.append([coord for coord in list(atom.get_coord())])
                last_atoms.append(coord[-1])

    args = []
    for i in range(l1[0], l1[-1] + 1):
        args.append(Point3D(last_atoms[i]))
    first_layer = Plane(*args)

    args = []
    for i in range(l2[0], l2[-1] + 1):
        args.append(Point3D(last_atoms[i]))
    second_layer = Plane(*args)

    args = []
    for i in range(l3[0], l3[-1] + 1):
        args.append(Point3D(last_atoms[i]))
    third_layer = Plane(*args)

    return [degrees(N(first_layer.angle_between(second_layer))), degrees(N(first_layer.angle_between(third_layer))),
            degrees(second_layer.angle_between(third_layer))]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='input_file',
                        required=True,
                        type=str)

    parser.add_argument('-l1', dest='fst_l',
                        nargs='+',
                        required=True,
                        type=int)

    parser.add_argument('-l2', dest='sd_l',
                        nargs='+',
                        required=True,
                        type=int)

    parser.add_argument('-l3', dest='td_l',
                        nargs='+',
                        required=True,
                        type=int)

    args = parser.parse_args()

    in_file_path = args.input_file
    first_lay = args.fst_l
    second_lay = args.sd_l
    third_lay = args.td_l

    plane_angle = plane_angle(in_file_path, first_lay, second_lay, third_lay)
    cols = ['prot_name', 'angle_between_first_second',
            'angle_between_first_third', 'angle_between_second_third']

    prot_name = path.basename(in_file_path)
    data = [prot_name]
    for elem in plane_angle:
        data.append(elem)

    df = pd.DataFrame([data], columns=cols)
    df.to_csv('planes_angles.csv', mode='a', header=False)
