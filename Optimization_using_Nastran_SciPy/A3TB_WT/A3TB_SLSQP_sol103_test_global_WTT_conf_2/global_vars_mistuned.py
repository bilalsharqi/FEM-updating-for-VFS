import numpy as np
import pandas as pd
import os
import shutil
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pyNastran.bdf.bdf import BDF
from pyMSCNastranUtilities import *

def init():
    global  mistuned_grids, mistuned_n_grids, mistuned_grid_coords, mistuned_freq_NASTRAN, mistuned_mode_shapes, mistuned_static_deform, M_mistuned, x_G_mistuned, J_G_mistuned

    # Read mistuned results
    file_path = "out"
    f06_file = 'mistuned_sol103.f06'
    model_coords = 'mistuned_sol103_coor.txt'
    n_modes=10
    n_subcases=1
    
    # read the reference result files and store the numerical data
    mistuned_grids, mistuned_n_grids, mistuned_grid_coords = importGrids(file_path, ['mistuned_A3TB1.bdf',model_coords], debug=True)
    mistuned_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
    mistuned_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, mistuned_n_grids, mistuned_grids, n_subcases,[],debug=True)
    mistuned_static_deform= importDisplacements(file_path, f06_file, n_subcases, mistuned_grids, grids_order=[], debug=True)
    M_mistuned, x_G_mistuned, J_G_mistuned = importRigidBodyMassData(file_path, f06_file,debug=True)
    print("\nMistuned NASTRAN data import completed")