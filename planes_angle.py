from Bio.PDB import *
from sympy import *
from os import path
from math import degrees
from COM_helix import COM_helix
import numpy as np
import argparse
import pandas as pd


def plane_angle(pdb_file, l1, l2, l3):
    com_hels = COM_helix(pdb_file)

    args = []
    for i in range(l1[0], l1[-1] + 1):
        args.append(Point3D(com_hels[i-1][0]))
    first_layer = Plane(*args)

    args = []
    for i in range(l2[0], l2[-1] + 1):
        args.append(Point3D(com_hels[i-1][0]))
    second_layer = Plane(*args)

    args = []
    for i in range(l3[0], l3[-1] + 1):
        args.append(Point3D(com_hels[i-1][0]))
    third_layer = Plane(*args)

    return [degrees(N(first_layer.angle_between(second_layer))), degrees(N(first_layer.angle_between(third_layer))),
            degrees(second_layer.angle_between(third_layer))]


print(plane_angle('/home/stephen/Desktop/PDB/1db1.pdb', [1, 2, 3], [4, 5, 6], [7, 8, 9]))
