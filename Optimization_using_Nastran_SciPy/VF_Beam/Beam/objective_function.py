# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 14:06:29 2021

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


# Objective function

def obj_func(x,inp_dir, run_dir, out_dir, ref_file,mistuned_file, M_ref, J_G_ref, x_G_ref, ref_freq_NASTRAN):
    
    # create an output directory for the results if one doesn't exist
    outputDir = "mistune_output"
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)
    
    # subroutine to update Nastran input for x (stiffness and mass properties)
    mistuned_model = BDF()
    # mistuned_model.read_bdf("mistune_output/mistunedBeam_Alt_1.bdf", punch=True)
    # print(mistuned_model.get_bdf_stats())
    mistuned_model.read_bdf("inp/mistunedBeam.bdf", punch=True)
    print(mistuned_model.get_bdf_stats())
       
    # get the number of properties of each type - separate mass and stiffness
    elemPropKeys = list(mistuned_model.properties.keys())
    # matStiffPropKeys = list(mistuned_model.materials.keys())
    # matMassPropKeys = list(mistuned_model.materials.keys())
    # conMassKeys = list(mistuned_model.masses.keys())
       
    # numElemProps = len(elemPropKeys)
    # numMatStiffProps = len(matStiffPropKeys)
    # numMatMassProps = len(matMassPropKeys)
    # numConMasses = len(conMassKeys)

    # Open BDF
    # model = BDF()
    # write material properties from updated design variables
    print(x)
    for i in range(len(elemPropKeys)):
        mistuned_model.materials[i+1].rho = x[i]
        mistuned_model.materials[i+1].e   = x[i+len(elemPropKeys)]
        # mistuned_model.materials[i+1].nu  = x[2][i]
        
    if mistuned_file.endswith('.bdf'):
        mistuned_file = mistuned_file[:-4]
    bdfFilenameOut = outputDir + '/' + mistuned_file + str(1) + '.bdf'
    mistuned_model.write_bdf(bdfFilenameOut,interspersed = False)
    
    
    # # Call Nastran for generating reference results
    # runNastran(inp_dir, run_dir, out_dir, ref_file, debug=True)
    
    # # Read reference results
    # # Note: There are 4 loading subcases and 15 eigenvalues computed for each
    # # deformed modal analysis
    # file_path = "out"
    # f06_file = 'sol400.f06'
    # model_coords = 'sol400_coor.txt'
    # n_modes=15
    # n_subcases=4
    
    # # Components
    
    # # read the reference result files and store the numerical data
    # ref_grids, ref_n_grids, ref_grid_coords = importGrids(file_path, ['refBeam.bdf',model_coords], debug=True)
    # ref_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
    # ref_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, ref_n_grids, ref_grids, n_subcases,[],debug=True)
    # ref_static_deform= importDisplacements(file_path, f06_file, n_subcases, ref_grids, grids_order=[], debug=True)
    # M_ref, x_G_ref, J_G_ref = importRigidBodyMassData(file_path, f06_file,debug=True)
    # print("\nReference NASTRAN data import completed")

    # Call Nastran for generating mistuned results
    runNastran("mistune_output", run_dir, out_dir, "mistuned_sol400.dat", debug=True)
    
    # Read mistuned results
    # Note: There are 4 loading subcases and 15 eigenvalues computed for each
    # deformed modal analysis
    file_path = "out"
    f06_file = 'mistuned_sol400.f06'
    model_coords = 'mistuned_sol400_coor.txt'
    n_modes=15
    n_subcases=4
    
    # read the reference result files and store the numerical data
    mistuned_grids, mistuned_n_grids, mistuned_grid_coords = importGrids(file_path, ['mistunedBeam1.bdf',model_coords], debug=True)
    mistuned_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
    mistuned_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, mistuned_n_grids, mistuned_grids, n_subcases,[],debug=True)
    mistuned_static_deform= importDisplacements(file_path, f06_file, n_subcases, mistuned_grids, grids_order=[], debug=True)
    M_mistuned, x_G_mistuned, J_G_mistuned = importRigidBodyMassData(file_path, f06_file,debug=True)
    print("\nMistuned NASTRAN data import completed")

    # Mass module
    # Load reference mass
    m_0 = M_ref
    
    # Load mistuned/current mass
    m = M_mistuned
    
    # Calculate mass module metric
    mass_metric = ((m-m_0)/m_0)**2
    
    # Inertia module 
    # Load reference inertia
    Ixx_0 = J_G_ref[0][0]
    Iyy_0 = J_G_ref[1][1]
    Izz_0 = J_G_ref[2][2]
    
    # Load mistuned/current inertia
    Ixx = J_G_mistuned[0][0]
    Iyy = J_G_mistuned[1][1]
    Izz = J_G_mistuned[2][2]
    
    # Calculate inertia module metric
    Ixx_metric = ((Ixx - Ixx_0)/Ixx_0)**2
    Iyy_metric = ((Iyy - Iyy_0)/Iyy_0)**2
    Izz_metric = ((Izz - Izz_0)/Izz_0)**2
    
    inertia_metric = Ixx_metric + Iyy_metric + Izz_metric
    
    # c.g. module
    # Load reference c.g.
    xcg_0 = x_G_ref[0]
    ycg_0 = x_G_ref[1]
    zcg_0 = x_G_ref[2]
    
    # Load mistuned/current c.g
    xcg = x_G_mistuned[0]
    ycg = x_G_mistuned[1]
    zcg = x_G_mistuned[2]
    
    # Calculate c.g. module metric
    # not dividing by the reference c.g. since 
    # it can be zero. BS - revisit this
    xcg_metric = ((xcg - xcg_0))**2
    ycg_metric = ((ycg - ycg_0))**2
    zcg_metric = ((zcg - zcg_0))**2
    
    cg_metric = xcg_metric + ycg_metric + zcg_metric
    
    # Frequency module
    # Load reference frequencies
    omega_0 = ref_freq_NASTRAN
    
    # Load mistuned frequencies
    omega = mistuned_freq_NASTRAN
    
    # Calculate frequency module metric
    # temp freq_metric variable to store individual frequency
    # metrics for all the frequencies for all the subcases
    freq_metric = [[np.zeros(1) for i in range(len(ref_freq_NASTRAN[0]))] for j in range(len(ref_freq_NASTRAN))]
    for i in range(len(ref_freq_NASTRAN)):
        for j in range(len(ref_freq_NASTRAN[0])):
            freq_metric[i][j] = ((omega[i][j]-omega_0[i][j])/omega_0[i][j])**2
    frequency_metric = np.sum(np.sum(freq_metric,0),0)
    
    obj = mass_metric + inertia_metric + cg_metric + frequency_metric
    
    return obj, mass_metric, inertia_metric, cg_metric, frequency_metric

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

# # one output objective function
# def objective_function(x):
#     outputs = obj_func(x,"inp", "run", "out", "sol400.dat","mistunedBeam_Alt_1.bdf", M_ref, J_G_ref, x_G_ref, ref_freq_NASTRAN)
#     objective = outputs[0]
#     return objective

# test_obj = obj_func(x,"inp", "run", "out", "sol400.dat","mistunedBeam.bdf")
# test_objective = objective_function(x)