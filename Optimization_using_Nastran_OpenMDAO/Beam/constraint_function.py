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

# plt.close("all")

from pyMSCNastranUtilities import *

# bnds = ((min_area,None),(min_area,None),(min_area,None),(min_area,None),(min_area,None),(min_area,None),(min_area,None),(min_area,None),(min_area,None),(min_area,None))
# Read reference results
# Note: There are 4 loading subcases and 15 eigenvalues computed for each
# deformed modal analysis
file_path = "out"
f06_file = 'sol400.f06'
model_coords = 'sol400_coor.txt'
n_modes=15
n_subcases=4

# read the reference result files and store the numerical data
ref_grids, ref_n_grids, ref_grid_coords = importGrids(file_path, ['refBeam.bdf',model_coords], debug=True)
ref_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
ref_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, ref_n_grids, ref_grids, n_subcases,[],debug=True)
ref_static_deform= importDisplacements(file_path, f06_file, n_subcases, ref_grids, grids_order=[], debug=True)
M_ref, x_G_ref, J_G_ref = importRigidBodyMassData(file_path, f06_file,debug=True)
print("\nReference NASTRAN data import completed")

# set up limits for the inequality constraints
ineq_const_limits = 1.1*np.float64([M_ref, J_G_ref[0][0], J_G_ref[1][1], J_G_ref[2][2],
                     x_G_ref[0],x_G_ref[1],x_G_ref[2]])

cons = (        {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] - constraint_func()[0]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func()[0]}, 
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[1] - constraint_func()[1]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func()[0]}, 
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[2] - constraint_func()[2]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func()[0]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[3] - constraint_func()[3]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func()[0]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[4] - constraint_func()[4]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func()[0]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[5] - constraint_func()[5]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func()[0]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[6] - constraint_func()[6]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func()[0]},
        )

# Constraint function

def constraint_func(x):
  
    # Read mistuned results
    # Note: There are 4 loading subcases and 15 eigenvalues computed for each
    # deformed modal analysis
    file_path = "out"
    f06_file = 'mistuned_sol400.f06'
    model_coords = 'mistuned_sol400_coor.txt'
    n_modes=15
    n_subcases=4
    
    # read the reference result files and store the numerical data
    mistuned_grids, mistuned_n_grids, mistuned_grid_coords = importGrids(file_path, ['mistunedBeam.bdf',model_coords], debug=True)
    mistuned_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
    mistuned_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, mistuned_n_grids, mistuned_grids, n_subcases,[],debug=True)
    mistuned_static_deform= importDisplacements(file_path, f06_file, n_subcases, mistuned_grids, grids_order=[], debug=True)
    M_mistuned, x_G_mistuned, J_G_mistuned = importRigidBodyMassData(file_path, f06_file,debug=True)
    print("\nMistuned NASTRAN data import completed")
    
    # Mass module
    # Load mistuned/current mass
    m = M_mistuned
    
    # Inertia module 
    # Load mistuned/current inertia
    Ixx = J_G_mistuned[0][0]
    Iyy = J_G_mistuned[1][1]
    Izz = J_G_mistuned[2][2]
    
    # c.g. module  
    # Load mistuned/current c.g
    xcg = x_G_mistuned[0]
    ycg = x_G_mistuned[1]
    zcg = x_G_mistuned[2]
    
    # MAC module
    MAC = [[np.zeros(1) for i in range(n_modes)] for j in range(n_subcases)]
    data_MAC_reference = [[np.zeros([6*mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]
    data_MAC_mistuned  = [[np.zeros([6*mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]
    for i in range(n_subcases):
        for j in range(n_modes):
            data_MAC_reference[i][j] = ref_mode_shapes[i][j].flatten('F')
            data_MAC_mistuned[i][j] = mistuned_mode_shapes[i][j].flatten('F')       
    
    test_MAC2 = ComputeMAC(data_MAC_reference[0][0],data_MAC_mistuned[0][0])
    for i in range(n_subcases):
        for j in range(n_modes):
                MAC[i][j] = ComputeMAC(data_MAC_reference[i][j], data_MAC_mistuned[i][j])
    
    return m,Ixx,Iyy,Izz,xcg,ycg,zcg,MAC

# Load the mistuned model bdf to provide initial guess
mistuned_model = BDF()
mistuned_model.read_bdf("inp/mistunedBeam.bdf", punch=True)
# print(mistuned_model.get_bdf_stats())
   
# get the number of properties of each type - separate mass and stiffness
elemPropKeys = list(mistuned_model.properties.keys())
matStiffPropKeys = list(mistuned_model.materials.keys())
matMassPropKeys = list(mistuned_model.materials.keys())
conMassKeys = list(mistuned_model.masses.keys())
   
numElemProps = len(elemPropKeys)
numMatStiffProps = len(matStiffPropKeys)
numMatMassProps = len(matMassPropKeys)
numConMasses = len(conMassKeys)

# read mistuned stiffness values
mistuned_stiffness = np.zeros(numMatStiffProps)
mistuned_mass = np.zeros(numMatMassProps)
mistuned_nu = np.zeros(numMatStiffProps)

for stiffness in range(len(matStiffPropKeys)):
    mistuned_stiffness[stiffness] = mistuned_model.materials[stiffness+1].e
    
# read mistuned mass values
for mass in range(len(matMassPropKeys)):
    mistuned_mass[mass] = mistuned_model.materials[mass+1].rho
    
# read reference nu values
for nu in range(len(matMassPropKeys)):
    mistuned_nu[nu] = mistuned_model.materials[nu+1].nu

x = np.array([np.transpose(mistuned_mass)])
x = np.append(x,[np.transpose(mistuned_stiffness)])
# x = np.append(x,[np.transpose(mistuned_nu)])
    
test_constraint = constraint_func(x)