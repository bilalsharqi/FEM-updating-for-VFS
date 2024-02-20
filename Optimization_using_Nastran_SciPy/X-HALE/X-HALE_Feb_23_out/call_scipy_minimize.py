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

# Call Nastran for generating reference results
runNastran("inp", "run", "out", "sol400.dat", debug=True)
runNastran("mistune_output", "run", "out", "mistuned_sol400.dat", debug=True)

# import and initialize global variables containing the mistuned info
import global_vars_mistuned
global_vars_mistuned.init()

# Read reference results
# Note: There are 1 loading subcases and 10 eigenvalues computed for each
# deformed modal analysis
file_path = "out"
f06_file = 'sol400.f06'
model_coords = 'sol400_coor.txt'
n_modes=15
n_subcases=1

# Components

# read the reference result files and store the numerical data
ref_grids, ref_n_grids, ref_grid_coords = importGrids(file_path, ['input_XHALE.bdf',model_coords], debug=True)
ref_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
ref_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, ref_n_grids, ref_grids, n_subcases,[],debug=True)
ref_static_deform= importDisplacements(file_path, f06_file, n_subcases, ref_grids, grids_order=[], debug=True)
M_ref, x_G_ref, J_G_ref = importRigidBodyMassData(file_path, f06_file,debug=True)
print("\nReference NASTRAN data import completed")
# Bilal - need to ensure the M_ref and others are passed correctly
# as arguments to the objective function so we do not run 
# the reference NASTRAN FEM every iteration

# force reference frequencies to match the experimentally identified ones
# Note: Not all modes are experimentally identified, so the other modes
# still come from the reference FEM
# values for subcase 1 - here it is the outboard config for X-HALE
ref_freq_NASTRAN[0][6] = 0.97
ref_freq_NASTRAN[0][7] = 2.73
ref_freq_NASTRAN[0][8] = 4.29
ref_freq_NASTRAN[0][9] = 4.53
ref_freq_NASTRAN[0][10] = 6.2
ref_freq_NASTRAN[0][11] = 7.22
ref_freq_NASTRAN[0][12] = 8.29

# values for subcase 2 - here it is the inboard config for X-HALE
# ref_freq_NASTRAN[0][6] = 0.95
# ref_freq_NASTRAN[0][7] = 2.7
# ref_freq_NASTRAN[0][8] = 4.5
# # ref_freq_NASTRAN[0][9] = 4.9
# # ref_freq_NASTRAN[0][10] = 5.9
# ref_freq_NASTRAN[0][12] = 7.3
# ref_freq_NASTRAN[0][13] = 8.3

# one output objective function
def objective_function(x):
    outputs = obj_func(x,"inp", "run", "out", "sol400.dat","mistuned_input_XHALE.bdf", M_ref, J_G_ref, x_G_ref, ref_freq_NASTRAN,x_lb,x_ub)
    objective = outputs[0]
    return objective
    
# Load the mistuned model bdf to provide initial guess
mistuned_model = BDF()
mistuned_model.read_bdf("inp/mistuned_input_XHALE.bdf", punch=True)
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
mistuned_stiffness = np.zeros(24)
# mistuned_mass = np.zeros(numMatMassProps)
mistuned_spring = np.zeros(6)
# mistuned_nu = np.zeros(numMatStiffProps)
# mistuned_thickness = np.zeros(numElemProps)
# mistuned_height = np.zeros(numElemProps)

# manually add initial stiffness variables to mistuned_stiffness
mistuned_stiffness[0] = mistuned_model.materials[5551].e
mistuned_stiffness[1] = mistuned_model.materials[5552].e
mistuned_stiffness[2] = mistuned_model.materials[5553].e
mistuned_stiffness[3] = mistuned_model.materials[5554].e
mistuned_stiffness[4] = mistuned_model.materials[5555].e
mistuned_stiffness[5] = mistuned_model.materials[5556].e
mistuned_stiffness[6] = mistuned_model.materials[6000].e
mistuned_stiffness[7] = mistuned_model.materials[6001].e
mistuned_stiffness[8] = mistuned_model.materials[6002].e
mistuned_stiffness[9] = mistuned_model.materials[6003].e
mistuned_stiffness[10] = mistuned_model.materials[6004].e
mistuned_stiffness[11] = mistuned_model.materials[1111].e

mistuned_stiffness[12] = mistuned_model.materials[5551].g
mistuned_stiffness[13] = mistuned_model.materials[5552].g
mistuned_stiffness[14] = mistuned_model.materials[5553].g
mistuned_stiffness[15] = mistuned_model.materials[5554].g
mistuned_stiffness[16] = mistuned_model.materials[5555].g
mistuned_stiffness[17] = mistuned_model.materials[5556].g
mistuned_stiffness[18] = mistuned_model.materials[6000].g
mistuned_stiffness[19] = mistuned_model.materials[6001].g
mistuned_stiffness[20] = mistuned_model.materials[6002].g
mistuned_stiffness[21] = mistuned_model.materials[6003].g
mistuned_stiffness[22] = mistuned_model.materials[6004].g
mistuned_stiffness[23] = mistuned_model.materials[1111].g

# manually add initial density variables to mistuned_mass
# mistuned_mass[0] = mistuned_model.materials[1111].rho
# mistuned_mass[1] = mistuned_model.materials[5555].rho
# mistuned_mass[2] = mistuned_model.materials[6000].rho
# mistuned_mass[3] = mistuned_model.materials[55555].rho

# manually add tail connector spring constants to mistuned_spring
mistuned_spring[0] = mistuned_model.properties[55555].Ki[0]
mistuned_spring[1] = mistuned_model.properties[55555].Ki[1]
mistuned_spring[2] = mistuned_model.properties[55555].Ki[2]
mistuned_spring[3] = mistuned_model.properties[55555].Ki[3]
mistuned_spring[4] = mistuned_model.properties[55555].Ki[4]
mistuned_spring[5] = mistuned_model.properties[55555].Ki[5]

# initial x
# x = np.array([np.transpose(mistuned_mass)])
x = np.transpose(mistuned_stiffness)
x = np.append(x,[np.transpose(mistuned_spring)])
# x = np.append(x,[np.transpose(mistuned_thickness)])
# x = np.append(x,[np.transpose(mistuned_height)])

# Load the reference model bdf using pynastran
# open the .bdf file
ref_model = BDF()
ref_model.read_bdf("inp/input_XHALE.bdf", punch=True)
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
ref_stiffness = np.zeros(24)
# ref_mass = np.zeros(numMatMassProps)
ref_spring = np.zeros(6)
# ref_nu = np.zeros(numMatStiffProps)
# ref_thickness = np.zeros(numElemProps)
# ref_height = np.zeros(numElemProps)

# manually add reference stiffness variables ref_stiffness
ref_stiffness[0] = ref_model.materials[5551].e
ref_stiffness[1] = ref_model.materials[5552].e
ref_stiffness[2] = ref_model.materials[5553].e
ref_stiffness[3] = ref_model.materials[5554].e
ref_stiffness[4] = ref_model.materials[5555].e
ref_stiffness[5] = ref_model.materials[5556].e
ref_stiffness[6] = ref_model.materials[6000].e
ref_stiffness[7] = ref_model.materials[6001].e
ref_stiffness[8] = ref_model.materials[6002].e
ref_stiffness[9] = ref_model.materials[6003].e
ref_stiffness[10] = ref_model.materials[6004].e
ref_stiffness[11] = ref_model.materials[1111].e

ref_stiffness[12] = ref_model.materials[5551].g
ref_stiffness[13] = ref_model.materials[5552].g
ref_stiffness[14] = ref_model.materials[5553].g
ref_stiffness[15] = ref_model.materials[5554].g
ref_stiffness[16] = ref_model.materials[5555].g
ref_stiffness[17] = ref_model.materials[5556].g
ref_stiffness[18] = ref_model.materials[6000].g
ref_stiffness[19] = ref_model.materials[6001].g
ref_stiffness[20] = ref_model.materials[6002].g
ref_stiffness[21] = ref_model.materials[6003].g
ref_stiffness[22] = ref_model.materials[6004].g
ref_stiffness[23] = ref_model.materials[1111].g


# manually add reference stiffness variables to ref_stiffness
# ref_mass[0] = ref_model.materials[1111].rho
# ref_mass[1] = ref_model.materials[5555].rho
# ref_mass[2] = ref_model.materials[6000].rho
# ref_mass[3] = ref_model.materials[55555].rho

# manually add root spring constants to ref_spring
ref_spring[0] = ref_model.properties[55555].Ki[0]
ref_spring[1] = ref_model.properties[55555].Ki[1]
ref_spring[2] = ref_model.properties[55555].Ki[2]
ref_spring[3] = ref_model.properties[55555].Ki[3]
ref_spring[4] = ref_model.properties[55555].Ki[4]
ref_spring[5] = ref_model.properties[55555].Ki[5]

# initial x from reference FEM
# x_ref = np.array([np.transpose(ref_mass)])
x_ref = np.transpose(ref_stiffness)
x_ref = np.append(x_ref,[np.transpose(ref_spring)])
# x_ref = np.append(x_ref,[np.transpose(ref_thickness)])
# x_ref = np.append(x_ref,[np.transpose(ref_height)])

# define bounds on design variables
# Bilal 2/6/23 these were commented out because we wanted to implement 
# separate bounds on each DV
# bnds_mass,x_lb_mass,x_ub_mass = bounds_func(x_ref[0:4],0.5,2.0)
bnds_stiffness,x_lb_stiffness,x_ub_stiffness = bounds_func(x_ref[0:24],0.75,1.25)
bnds_spring,x_lb_spring,x_ub_spring = bounds_func(x_ref[24:],0.5,2.0)

# separate bounds on each design variable
# bnds_mass_1,x_lb_mass_1,x_ub_mass_1 = bounds_func_single(x_ref[0],0.9,1.1)
# bnds_mass_2,x_lb_mass_2,x_ub_mass_2 = bounds_func_single(x_ref[1],0.9,1.1)
# bnds_mass_3,x_lb_mass_3,x_ub_mass_3 = bounds_func_single(x_ref[2],0.9,1.1)
# bnds_stiffness_1,x_lb_stiffness_1,x_ub_stiffness_1 = bounds_func_single(x_ref[3],0.75,1.25)
# bnds_stiffness_2,x_lb_stiffness_2,x_ub_stiffness_2 = bounds_func_single(x_ref[4],0.5,5.0)
# bnds_stiffness_3,x_lb_stiffness_3,x_ub_stiffness_3 = bounds_func_single(x_ref[5],0.95,1.05)
# bnds_stiffness_4,x_lb_stiffness_4,x_ub_stiffness_4 = bounds_func_single(x_ref[6],0.95,1.05)
# bnds_spring_1,x_lb_spring_1,x_ub_spring_1 = bounds_func_single(x_ref[7],0.5,2.0)
# bnds_spring_2,x_lb_spring_2,x_ub_spring_2 = bounds_func_single(x_ref[8],0.5,2.0)
# bnds_spring_3,x_lb_spring_3,x_ub_spring_3 = bounds_func_single(x_ref[9],0.0,2.0)
# bnds_spring_4,x_lb_spring_4,x_ub_spring_4 = bounds_func_single(x_ref[10],0.5,2.0)
# bnds_spring_5,x_lb_spring_5,x_ub_spring_5 = bounds_func_single(x_ref[11],0.5,2.0)

# assemble lower bounds 
# x_lb_mass = np.array([x_lb_mass_1, x_lb_mass_2, x_lb_mass_3])
# x_lb_stiffness = np.array([x_lb_stiffness_1, x_lb_stiffness_2, x_lb_stiffness_3, x_lb_stiffness_4])
# x_lb_spring = [x_lb_spring_1, x_lb_spring_2, x_lb_spring_3, x_lb_spring_4, x_lb_spring_5]

# assemble upper bounds
# x_ub_mass = np.array([x_ub_mass_1, x_ub_mass_2, x_ub_mass_3])
# x_ub_stiffness = np.array([x_ub_stiffness_1, x_ub_stiffness_2, x_ub_stiffness_3, x_ub_stiffness_4])
# x_ub_spring = [x_ub_spring_1, x_ub_spring_2, x_ub_spring_3, x_ub_spring_4, x_ub_spring_5]

x_lb = np.array(x_lb_stiffness)
# x_lb = np.append(x_lb,x_lb_stiffness)
x_lb = np.append(x_lb,x_lb_spring)

x_ub = np.array(x_ub_stiffness)
# x_ub = np.append(x_ub,x_ub_stiffness)
x_ub = np.append(x_ub,x_ub_spring)

# scaling design variables to be between 0 and 1 each
x_new = (x - x_lb)/(x_ub - x_lb)

# creating new bounds for DV
bnds = bounds_func(np.ones(30),0.0,1.0)[0]

# set up limits for the inequality constraints
ineq_const_limits_lb = 0.9*np.float64([M_ref, J_G_ref[0][0], J_G_ref[1][1], J_G_ref[2][2],
                      x_G_ref[0],x_G_ref[1],x_G_ref[2]])

ineq_const_limits_ub = 1.1*np.float64([M_ref, J_G_ref[0][0], J_G_ref[1][1], J_G_ref[2][2],
                      x_G_ref[0],x_G_ref[1],x_G_ref[2]])

# append the inequality constraint for MAC
ineq_const_limits_lb = np.append(ineq_const_limits_lb, 0.8)
ineq_const_limits_ub = np.append(ineq_const_limits_ub, 1.0)


cons = (        {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[0] - constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[0]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[0] + constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[0]}, 
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[1] - constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[1]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[1] + constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[1]}, 
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[2] - constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[2]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[2] + constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[2]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[3] - constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[3]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[3] + constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[3]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[4] - constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[4]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[4] + constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[4]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[5] - constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[5]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[5] + constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[5]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[6] - constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[6]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[6] + constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[6]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits_ub[7] - constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[7][7]}, {'type': 'ineq', 'fun': lambda x: -ineq_const_limits_lb[7] + constraint_func(x, ref_grids, ref_n_grids, ref_grid_coords, ref_freq_NASTRAN, ref_mode_shapes, ref_static_deform, M_ref, x_G_ref, J_G_ref)[7][7]},
        )
 
# unconstrained Nelder-Mead optimization
scipy_uncon_nm = sio.minimize(objective_function, x_new, args=(),method='Nelder-Mead', bounds=bnds,options={'return_all': True, 'fatol': 1e-2,'disp': True, 'adaptive': True})        
print("SciPy unconstrained NM:", scipy_uncon_nm)
x_opt_uncon = scipy_uncon_nm.x

# constrained SLSQP optimization
# scipy_con_sqp  = sio.minimize(objective_function, x_new, method='SLSQP', jac='cs', bounds=bnds, constraints=cons, tol=1e-6, options={'maxiter': 500,'ftol': 1e-6,'disp':True, 'eps': 3e-2})
# print("SciPy constrained SQP:", scipy_con_sqp)
# x_opt_con   = scipy_con_sqp.x
