# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 12:00:26 2021

@author: bilal
"""


import numpy as np
import os
import scipy.optimize as sio
import shutil
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib import cm
from pyNastran.bdf.bdf import BDF

# plt.close("all")

from pyMSCNastranUtilities import *
from objective_function import *
from constraint_function import *
from bounds_function import *

# using SciPy minimize 

# Call Nastran for generating reference results
# runNastran("inp", "run", "out", "sol103.dat", debug=True)
runNastran("mistune_output", "run", "out", "mistuned_sol103.dat", debug=True)
# Read reference results
# Note: There are 4 loading subcases and 15 eigenvalues computed for each
# deformed modal analysis
file_path = "out"
f06_file = 'sol103.f06'
model_coords = 'sol103_coor.txt'
n_modes=15
n_subcases=1

# Components

# read the reference result files and store the numerical data
ref_grids, ref_n_grids, ref_grid_coords = importGrids(file_path, ['refBeam.bdf',model_coords], debug=True)
ref_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
ref_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, ref_n_grids, ref_grids, n_subcases,[],debug=True)
ref_static_deform= importDisplacements(file_path, f06_file, n_subcases, ref_grids, grids_order=[], debug=True)
M_ref, x_G_ref, J_G_ref = importRigidBodyMassData(file_path, f06_file,debug=True)
print("\nReference NASTRAN data import completed")
# Bilal - need to ensure the M_ref and others are passed correctly
# as arguments to the objective function so we do not run 
# the reference NASTRAN FEM every iteration

# one output objective function
def objective_function(x):
    outputs = obj_func(x,"inp", "run", "out", "sol103.dat","mistunedBeam.bdf", M_ref, J_G_ref, x_G_ref, ref_freq_NASTRAN,x_lb,x_ub)
    objective = outputs[0]
    return objective
    
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

# initialize vectors for reading model properties
mistuned_stiffness = np.zeros(numMatStiffProps)
mistuned_mass = np.zeros(numMatMassProps)
mistuned_nu = np.zeros(numMatStiffProps)
mistuned_thickness = np.zeros(numElemProps)
mistuned_height = np.zeros(numElemProps)

# read initial stiffness values
for stiffness in range(len(matStiffPropKeys)):
    mistuned_stiffness[stiffness] = mistuned_model.materials[stiffness+1].e
    
# read initial mass values
for mass in range(len(matMassPropKeys)):
    mistuned_mass[mass] = mistuned_model.materials[mass+1].rho
    
# read initial nu values
for nu in range(len(matMassPropKeys)):
    mistuned_nu[nu] = mistuned_model.materials[nu+1].nu

# read initial thickness values
for thickness in range(len(elemPropKeys)):
    mistuned_thickness[thickness] = mistuned_model.properties[thickness+1].dim[0][0]

# read initial height values
for height in range(len(elemPropKeys)):
    mistuned_height[height] = mistuned_model.properties[height+1].dim[0][1]

# initial x
x = np.array([np.transpose(mistuned_mass)])
x = np.append(x,[np.transpose(mistuned_stiffness)])
x = np.append(x,[np.transpose(mistuned_thickness)])
x = np.append(x,[np.transpose(mistuned_height)])


# Load the reference model bdf using pynastran
# open the .bdf file
ref_model = BDF()
ref_model.read_bdf("inp/refBeam.bdf", punch=True)
# print(ref_model.get_bdf_stats())
   
# get the number of properties of each type - separate mass and stiffness
elemPropKeys = list(ref_model.properties.keys())
matStiffPropKeys = list(ref_model.materials.keys())
matMassPropKeys = list(ref_model.materials.keys())
conMassKeys = list(ref_model.masses.keys())
   
numElemProps = len(elemPropKeys)
numMatStiffProps = len(matStiffPropKeys)
numMatMassProps = len(matMassPropKeys)
numConMasses = len(conMassKeys)

# initialize vectors for reading model properties
ref_stiffness = np.zeros(numMatStiffProps)
ref_mass = np.zeros(numMatMassProps)
ref_nu = np.zeros(numMatStiffProps)
ref_thickness = np.zeros(numElemProps)
ref_height = np.zeros(numElemProps)

# read initial stiffness values
for stiffness in range(len(matStiffPropKeys)):
    ref_stiffness[stiffness] = ref_model.materials[stiffness+1].e
    
# read initial mass values
for mass in range(len(matMassPropKeys)):
    ref_mass[mass] = ref_model.materials[mass+1].rho
    
# read initial nu values
for nu in range(len(matMassPropKeys)):
    ref_nu[nu] = ref_model.materials[nu+1].nu

# read initial thickness values
for thickness in range(len(elemPropKeys)):
    ref_thickness[thickness] = ref_model.properties[thickness+1].dim[0][0]

# read initial height values
for height in range(len(elemPropKeys)):
    ref_height[height] = ref_model.properties[height+1].dim[0][1]

# initial x from reference FEM
x_ref = np.array([np.transpose(ref_mass)])
x_ref = np.append(x_ref,[np.transpose(ref_stiffness)])
x_ref = np.append(x_ref,[np.transpose(ref_thickness)])
x_ref = np.append(x_ref,[np.transpose(ref_height)])
    
# x_new = np.float64([1,1,1,1,1,1,1,1,1,1,1,1])

bnds_temp,x_lb,x_ub = bounds_func(x_ref,0.9,1.1)

# scaling design variables to be between 0 and 1 each
x_new = (x - x_lb)/(x_ub - x_lb)

# creating new bounds for DV
bnds = bounds_func(np.ones(12),0.0,1.0)[0]

# set up limits for the inequality constraints
ineq_const_limits_lb = 0.9*np.float64([M_ref, J_G_ref[0][0], J_G_ref[1][1], J_G_ref[2][2],
                      x_G_ref[0],x_G_ref[1],x_G_ref[2]])
ineq_const_limits_ub = 1.1*np.float64([M_ref, J_G_ref[0][0], J_G_ref[1][1], J_G_ref[2][2],
                      x_G_ref[0],x_G_ref[1],x_G_ref[2]])


cons = (        {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[0] - constraint_func(x)[0]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[0] + constraint_func(x)[0]}, 
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[1] - constraint_func(x)[1]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[1] + constraint_func(x)[1]}, 
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[2] - constraint_func(x)[2]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[2] + constraint_func(x)[2]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[3] - constraint_func(x)[3]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[3] + constraint_func(x)[3]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[4] - constraint_func(x)[4]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[4] + constraint_func(x)[4]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[5] - constraint_func(x)[5]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[5] + constraint_func(x)[5]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[6] - constraint_func(x)[6]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[6] + constraint_func(x)[6]},
        )
 
# scipy_uncon_nm = sio.minimize(objective_function, x_new, args=(),method='Nelder-Mead', bounds=bnds,options={'return_all': True, 'fatol': 1e-2,'disp': True, 'adaptive': True})        
# print("SciPy unconstrained NM:", scipy_uncon_nm)
# x_opt_uncon = scipy_uncon_nm.x

scipy_con_sqp  = sio.minimize(objective_function, x_new, method='SLSQP', jac=None, bounds=bnds, constraints=cons, tol=1e-6, options={'maxiter': 500,'ftol': 1e-6,'disp':True, 'eps': 1e-2})
print("SciPy constrained SQP:", scipy_con_sqp)
x_opt_con   = scipy_con_sqp.x
