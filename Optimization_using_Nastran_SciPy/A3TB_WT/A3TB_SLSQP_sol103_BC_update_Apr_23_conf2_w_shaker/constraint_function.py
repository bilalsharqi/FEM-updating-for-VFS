# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 22:31:46 2021

@author: bilal
"""

import numpy as np
import pandas as pd
import os
import shutil
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pyNastran.bdf.bdf import BDF
from pyMSCNastranUtilities import *
import global_vars_mistuned

# Constraint function

def constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref):
    
    # # Read reference results
    # file_path = "out"
    # f06_file = 'sol103.f06'
    # model_coords = 'sol103_coor.txt'
    # n_modes=10
    # n_subcases=1
    
    # # Components
    
    # # read the reference result files and store the numerical data
    # ref_grids, ref_n_grids, ref_grid_coords = importGrids(file_path, ['ref_A3TB.bdf',model_coords], debug=True)
    # ref_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
    # ref_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, ref_n_grids, ref_grids, n_subcases,[],debug=True)
    # ref_static_deform= importDisplacements(file_path, f06_file, n_subcases, ref_grids, grids_order=[], debug=True)
    # M_ref, x_G_ref, J_G_ref = importRigidBodyMassData(file_path, f06_file,debug=True)
    # print("\nReference NASTRAN data import completed")

    # # Read mistuned results
    # file_path = "out"
    # f06_file = 'mistuned_sol103.f06'
    # model_coords = 'mistuned_sol103_coor.txt'
    n_modes=10
    n_subcases=1
    
    # # read the reference result files and store the numerical data
    # mistuned_grids, mistuned_n_grids, mistuned_grid_coords = importGrids(file_path, ['mistuned_A3TB1.bdf',model_coords], debug=True)
    # mistuned_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
    # mistuned_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, mistuned_n_grids, mistuned_grids, n_subcases,[],debug=True)
    # mistuned_static_deform= importDisplacements(file_path, f06_file, n_subcases, mistuned_grids, grids_order=[], debug=True)
    # M_mistuned, x_G_mistuned, J_G_mistuned = importRigidBodyMassData(file_path, f06_file,debug=True)
    # print("\nMistuned NASTRAN data import completed")
    
    # Mass module
    # Load mistuned/current mass
    m = global_vars_mistuned.M_mistuned
    
    # Inertia module 
    # Load mistuned/current inertia
    Ixx = global_vars_mistuned.J_G_mistuned[0][0]
    Iyy = global_vars_mistuned.J_G_mistuned[1][1]
    Izz = global_vars_mistuned.J_G_mistuned[2][2]
    
    # c.g. module  
    # Load mistuned/current c.g
    xcg = np.float64(global_vars_mistuned.x_G_mistuned[0])
    ycg = np.float64(global_vars_mistuned.x_G_mistuned[1])
    zcg = np.float64(global_vars_mistuned.x_G_mistuned[2])
    
    # MAC module
    temp_MAC_size = (n_subcases,n_subcases)
    
    MAC = [np.zeros(temp_MAC_size) for j in range(n_subcases)]
    
    data_MAC_reference = [[np.zeros([6*global_vars_mistuned.mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]
    data_MAC_mistuned  = [[np.zeros([6*global_vars_mistuned.mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]
    
    for i in range(n_subcases):
        for j in range(n_modes):
            data_MAC_reference[i][j] = ref_mode_shapes[i][j].flatten('C')
            data_MAC_mistuned[i][j] = global_vars_mistuned.mistuned_mode_shapes[i][j].flatten('C')
    
    for i in range(n_subcases):
            MAC[i] = mac(np.array(data_MAC_reference[i][:]), np.array(data_MAC_mistuned[i][:]))
    
    return m,Ixx,Iyy,Izz,xcg,ycg,zcg,MAC
    
