# generate_contact_matrices.py
#
# Usage:
#     python generate_contact_matrices.py input_file.csv output_file.mat
#
# Reads an edge list file, and generates a set of contact matrices in a Matlab
# data file.  The matlab file will contain (currently), 4 matrices:
#     C - couse contact matrix
#     M - Movement between classes contact matrix
#     R - Residence contact matrix
#     S - Residence section contact matrix
# These matrices are all counts of events of that type, and should be scaled
# as appropriate to represent numbers of epidemiological contacts.

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.io as sio
from covid19_Net import covid19_Net
import graph_ops
from covid_simulation import repeated_runs
import datetime
import visualization as viz
import sys
import pandas as pd

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    data = graph_ops.loadEdgeList(input_file)
    C,M,R,S = graph_ops.edgesAsMatrices(data)

    matrix_data = {
        'C' : C,
        'M' : M,
        'R' : R,
        'S' : S
        }

    sio.savemat(output_file, matrix_data)

if __name__== "__main__":
    main()

