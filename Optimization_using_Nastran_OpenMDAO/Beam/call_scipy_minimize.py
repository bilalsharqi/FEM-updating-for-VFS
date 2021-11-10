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

# initial x
x = np.array([np.transpose(mistuned_mass)])
x = np.append(x,[np.transpose(mistuned_stiffness)])
# x = np.append(x,[np.transpose(mistuned_nu)])

# bnds = ((2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),(2430,2970),
#         (6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700),(6300,7700))

bnds = bounds_func(0.9,1.1)

cons = (        {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] - constraint_func(x)[0]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func(x)[0]}, 
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[1] - constraint_func(x)[1]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func(x)[0]}, 
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[2] - constraint_func(x)[2]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func(x)[0]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[3] - constraint_func(x)[3]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func(x)[0]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[4] - constraint_func(x)[4]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func(x)[0]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[5] - constraint_func(x)[5]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func(x)[0]},
                {'type': 'ineq', 'fun': lambda x: ineq_const_limits[6] - constraint_func(x)[6]}, {'type': 'ineq', 'fun': lambda x: ineq_const_limits[0] + constraint_func(x)[0]},
        )
 
scipy_uncon_nm = sio.minimize(objective_function, x, args=(),method='Nelder-Mead', bounds=bnds,options={'xatol': 0.0001, 'fatol': 0.0001,'gtol': 1e-6, 'disp': True})        
# scipy_con_sqp  = sio.minimize(objective_function, x, method='SLSQP', jac=None, bounds=bnds, constraints=cons, tol=1e-6, options={'maxiter': 500,'ftol': 1e-6,'disp':True})

print("SciPy unconstrained NM:", scipy_uncon_nm)
# print("SciPy constrained SQP:", scipy_con_sqp)

x_opt_uncon = scipy_uncon_nm.x
# x_opt_con   = scipy_con_sqp.x
