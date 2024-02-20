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

from pyMSCNastranUtilities import *
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
    
    # Loop over design variables
    for i in range(n_x):
        
        # relative step size
        delta_x = np.maximum(h,np.multiply(h,x[i]))
        
        # update the DV
        x[i] += delta_x
        
        # Loop over constraint variables
        for j in range(n_c): 

            # Module to write updated model bdf with (x+delta_x)

            # create an output directory for the results if one doesn't exist
            outputDir = "mistune_output"
            if not os.path.isdir(outputDir):
                os.mkdir(outputDir)
            
            # subroutine to update Nastran input for perturbed x
            # read the mistuned model and list/record properties
            mistuned_model = BDF()
            mistuned_model.read_bdf("inp/mistunedBeam.bdf", punch=True)
            
            # BS - temporarily hardcoding mistuned_file name variable here
            mistuned_file = "mistunedBeam.bdf"
            # print(mistuned_model.get_bdf_stats())
               
            # get the number of properties of each type - separate mass and stiffness
            # Note: BS - this is setup assuming the number of MAT and PROP
            # cards in the Nastran model are the same, so reading only one is fine!
            elemPropKeys = list(mistuned_model.properties.keys())

            # write material properties from updated design variables
            print(x)
            for ii in range(len(elemPropKeys)):
                mistuned_model.materials[ii+1].rho = x[ii]
                mistuned_model.materials[ii+1].e   = x[ii+len(elemPropKeys)]
                mistuned_model.properties[ii+1].dim[0][0] = x[ii+2*len(elemPropKeys)]
                mistuned_model.properties[ii+1].dim[0][1] = x[ii+3*len(elemPropKeys)]
                
            if mistuned_file.endswith('.bdf'):
                mistuned_file = mistuned_file[:-4]
            bdfFilenameOut = outputDir + '/' + mistuned_file + str(1) + '.bdf'
            mistuned_model.write_bdf(bdfFilenameOut,interspersed = False)
            
            # Call to run Nastran with updated x (+ delta_x)
            runNastran("mistune_output", "run", "out", "mistuned_sol400.dat", debug=True)
            
            # Read perturbed, mistuned results
            # Note: There are 4 loading subcases and 15 eigenvalues computed for each
            # deformed modal analysis
            file_path = "out"
            f06_file = 'mistuned_sol400.f06'
            model_coords = 'mistuned_sol400_coor.txt'
            n_modes=15
            n_subcases=4
            
            # read the result files and store the numerical data
            mistuned_grids, mistuned_n_grids, mistuned_grid_coords = importGrids(file_path, ['mistunedBeam1.bdf',model_coords], debug=True)
            mistuned_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
            mistuned_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, mistuned_n_grids, mistuned_grids, n_subcases,[],debug=True)
            mistuned_static_deform= importDisplacements(file_path, f06_file, n_subcases, mistuned_grids, grids_order=[], debug=True)
            M_mistuned, x_G_mistuned, J_G_mistuned = importRigidBodyMassData(file_path, f06_file,debug=True)
            print("\nMistuned NASTRAN data import completed")

            # Once results from each DV perturbation are available for 
            # each constraint variable, only then
            # does it make sense to compute derivatives!

            # Compute FD derivatives of each constraint wrt each DV
            constr_plus = np.float64(constraint_func(x)[0:7])
            grad_constr[i][j] = (constr_plus[j] - constr[j])/delta_x  
            
            # reset x
            x[i] -= delta_x    
    return grad_constr

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

# # initialize vectors for reading model properties
# mistuned_stiffness = np.zeros(numMatStiffProps)
# mistuned_mass = np.zeros(numMatMassProps)
# mistuned_nu = np.zeros(numMatStiffProps)
# mistuned_thickness = np.zeros(numElemProps)
# mistuned_height = np.zeros(numElemProps)

# # read initial stiffness values
# for stiffness in range(len(matStiffPropKeys)):
#     mistuned_stiffness[stiffness] = mistuned_model.materials[stiffness+1].e
    
# # read initial mass values
# for mass in range(len(matMassPropKeys)):
#     mistuned_mass[mass] = mistuned_model.materials[mass+1].rho
    
# # read initial nu values
# for nu in range(len(matMassPropKeys)):
#     mistuned_nu[nu] = mistuned_model.materials[nu+1].nu

# # read initial thickness values
# for thickness in range(len(elemPropKeys)):
#     mistuned_thickness[thickness] = mistuned_model.properties[thickness+1].dim[0][0]

# # read initial height values
# for height in range(len(elemPropKeys)):
#     mistuned_height[height] = mistuned_model.properties[height+1].dim[0][1]

# # initial x
# x = np.array([np.transpose(mistuned_mass)])
# x = np.append(x,[np.transpose(mistuned_stiffness)])
# x = np.append(x,[np.transpose(mistuned_thickness)])
# x = np.append(x,[np.transpose(mistuned_height)])
# # x = np.append(x,[np.transpose(mistuned_nu)])

# test_grad_constr = grad_constraint_fd(x)