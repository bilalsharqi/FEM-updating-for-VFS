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
    
    print("\nEntering constraint function")

    n_modes=15
    n_subcases=4
    
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
    
    print("\nLeaving constraint function")
    return m,Ixx,Iyy,Izz,xcg,ycg,zcg,MAC
    
