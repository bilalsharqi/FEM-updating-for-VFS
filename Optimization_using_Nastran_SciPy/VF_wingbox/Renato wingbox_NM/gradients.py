# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 11:46:22 2021

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

# from pyMSCNastranUtilities import *
# from objective_function import obj_func
from constraint_function import constraint_func


# # Gradient of objective function using forward finite difference
# def grad_obj_fd(x):         
#     f     = objective_function(x)
#     n_x   = len(x)
#     grad_f = np.array([np.zeros(n_x)])
#     h      = 10**(-6)
#     for i in range(n_x):
#         delta_x = np.maximum(h,np.multiply(h,x[i]))
#         x[i] += delta_x
#         f_plus = objective_function(x)
#         grad_f[i] = (f_plus - f)/delta_x  
#         x[i] -= delta_x    
#     return grad_f


# Gradient of constraint function using forward finite difference
def grad_constraint_fd(x):          
    constr   = np.float64(constraint_func(x)[0:7])
    n_x      = len(x)
    n_c      = len(constr)
    grad_constr = (np.zeros((n_x,n_c)))
    h         = 10**(-6)
    
    for i in range(n_x):
        
        for j in range(n_c): 
            delta_x = np.maximum(h,np.multiply(h,x[i]))
            x[i] += delta_x
            # Need to add module to write updated model bdf with (x+delta_x)
            # Also need to add a call to run Nastran with updated x
            # Once results from each DV perturbation are available for 
            # each constraint variable, only then
            # does it make sense to compute derivatives!
            constr_plus = np.float64(constraint_func(x)[0:7])
            grad_constr[i][j] = (constr_plus[j] - constr[j])/delta_x  
            x[i] -= delta_x    
    return grad_constr

# # Load the mistuned model bdf to provide initial guess
# mistuned_model = BDF()
# mistuned_model.read_bdf("inp/mistunedBeam.bdf", punch=True)
# # print(mistuned_model.get_bdf_stats())
   
# # get the number of properties of each type - separate mass and stiffness
# elemPropKeys = list(mistuned_model.properties.keys())
# matStiffPropKeys = list(mistuned_model.materials.keys())
# matMassPropKeys = list(mistuned_model.materials.keys())
# conMassKeys = list(mistuned_model.masses.keys())
   
# numElemProps = len(elemPropKeys)
# numMatStiffProps = len(matStiffPropKeys)
# numMatMassProps = len(matMassPropKeys)
# numConMasses = len(conMassKeys)

# # read mistuned stiffness values
# mistuned_stiffness = np.zeros(numMatStiffProps)
# mistuned_mass = np.zeros(numMatMassProps)
# mistuned_nu = np.zeros(numMatStiffProps)

# for stiffness in range(len(matStiffPropKeys)):
#     mistuned_stiffness[stiffness] = mistuned_model.materials[stiffness+1].e
    
# # read mistuned mass values
# for mass in range(len(matMassPropKeys)):
#     mistuned_mass[mass] = mistuned_model.materials[mass+1].rho
    
# # read reference nu values
# for nu in range(len(matMassPropKeys)):
#     mistuned_nu[nu] = mistuned_model.materials[nu+1].nu

# x = np.array([np.transpose(mistuned_mass)])
# x = np.append(x,[np.transpose(mistuned_stiffness)])
# # x = np.append(x,[np.transpose(mistuned_nu)])

# test_grad_constr = grad_constraint_fd(x)