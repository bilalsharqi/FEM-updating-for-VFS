# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 18:12:51 2021

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

plt.close("all")

from pyMSCNastranUtilities import *

# Script to perform I/O operations and call Nastran

# Prepare arguments to call Nastran for reference case

# paths for folders
inp_dir = "inp"
out_dir = "out"
run_dir = "run"

# Call Nastran for generating reference results
runNastran(inp_dir, run_dir, out_dir, "sol400.dat", debug=True)

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
    
# Plot the static deformation for the reference case
# fig=plt.figure(figsize=(15,8))
# plt.plot(ref_grid_coords[0,:] + ref_static_deform[0][0][:]-1, ref_static_deform[0][2][:],'r',label='Gravity load' )
# plt.plot(ref_grid_coords[0,:] + ref_static_deform[1][0][:]-1, ref_static_deform[1][2][:],'b',label='10 N distributed load')
# plt.plot(ref_grid_coords[0,:] + ref_static_deform[2][0][:]-1, ref_static_deform[2][2][:],'g',label='25 N distributed load')
# plt.plot(ref_grid_coords[0,:]-1, ref_grid_coords[2,:]-1, 'k',label='No load')
# plt.rcParams["lines.linewidth"] = 3
# plt.title('Reference Static deformation with loading',fontsize=32 )
# plt.xlabel("Beam length [m]",fontsize=26)
# plt.ylabel("Vertical displacement [m]",fontsize=26)
# plt.ax = plt.gca()
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# plt.legend(fontsize=22)
# # plt.ylim(-1,1.5)

# Save figure
# fig.savefig("out/beam_static_1_5m_z.svg",bbox_inches='tight')

# Call Nastran for generating mistuned results
runNastran(inp_dir, run_dir, out_dir, "mistuned_sol400.dat", debug=True)

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
   
# Plot the static deformation for the mistuned case 
# fig=plt.figure(figsize=(15,8))
# plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[0][0][:]-1, mistuned_static_deform[0][2][:],'r',label='Gravity load' )
# plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[1][0][:]-1, mistuned_static_deform[1][2][:],'b',label='10 N distributed load')
# plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[2][0][:]-1, mistuned_static_deform[2][2][:],'g',label='25 N distributed load')
# plt.plot(mistuned_grid_coords[0,:]-1, mistuned_grid_coords[2,:]-1, 'k',label='No load')
# plt.rcParams["lines.linewidth"] = 3
# plt.title('Mistuned Static deformation with loading',fontsize=32 )
# plt.xlabel("Beam length [m]",fontsize=26)
# plt.ylabel("Vertical displacement [m]",fontsize=26)
# plt.ax = plt.gca()
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# plt.legend(fontsize=22)
# plt.ylim(-1,1.5)

# Save figure
# fig.savefig("out/mistuned_beam_static_1_5m_z.svg",bbox_inches='tight')

# MAC computation between reference and mistuned cases
# initialize MAC ndarray for n_subcases*n_modes*n_dof 
def mac(Psi_1, Psi_2):
    """Modal Assurance Criterion.
    Parameters
    ----------
    Psi_1, Psi_2 : float arrays
        Mode shape matrices to be compared.
    Returns
    -------
    mac : float array
    Examples
    --------
    """
    nummodes = Psi_1.shape[1]
    MAC = np.zeros((nummodes, nummodes))
    if Psi_1.shape == Psi_2.shape:
        for i in np.arange(nummodes):
            for j in np.arange(nummodes):
                MAC[i, j] = (abs(np.conj(Psi_1[:, i]) @ Psi_2[:, j])**2 /
                             (np.conj(Psi_1[:, i]) @ Psi_1[:, i] *
                              np.conj(Psi_2[:, j]) @ Psi_2[:, j]))
    else:
        print('Mode shape arrays must have the same size.')
    return MAC

n_dof = 6
MAC = [[np.zeros(1) for i in range(n_modes)] for j in range(n_subcases)]
data_MAC_reference = [[np.zeros([n_dof*mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]
data_MAC_mistuned  = [[np.zeros([n_dof*mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]
for i in range(n_subcases):
    for j in range(n_modes):
        data_MAC_reference[i][j] = ref_mode_shapes[i][j].flatten('F')
        data_MAC_mistuned[i][j] = mistuned_mode_shapes[i][j].flatten('F')       

# test_MAC = mac(np.asarray(data_MAC_reference[0]),np.asarray(data_MAC_mistuned[0]))
test_MAC2 = ComputeMAC(data_MAC_reference[0][0],data_MAC_mistuned[0][0])
for i in range(n_subcases):
    for j in range(n_modes):
            MAC[i][j] = ComputeMAC(data_MAC_reference[i][j], data_MAC_mistuned[i][j])

# Prepare arguments to be sent to optimizer
# Load the reference model bdf using pynastran
# open the .bdf file
ref_model = BDF()
ref_model.read_bdf("inp/refBeam.bdf", punch=True)
print(ref_model.get_bdf_stats())
   
# get the number of properties of each type - separate mass and stiffness
elemPropKeys = list(ref_model.properties.keys())
matStiffPropKeys = list(ref_model.materials.keys())
matMassPropKeys = list(ref_model.materials.keys())
conMassKeys = list(ref_model.masses.keys())
   
numElemProps = len(elemPropKeys)
numMatStiffProps = len(matStiffPropKeys)
numMatMassProps = len(matMassPropKeys)
numConMasses = len(conMassKeys)

# read reference stiffness values
ref_stiffness = np.zeros(numMatStiffProps)
ref_mass = np.zeros(numMatMassProps)
ref_nu = np.zeros(numMatStiffProps)

for stiffness in range(len(matStiffPropKeys)-1):
    ref_stiffness[stiffness] = ref_model.materials[stiffness+1].e
    
# read reference mass values
for mass in range(len(matMassPropKeys)-1):
    ref_mass[mass] = ref_model.materials[mass+1].rho
    
# read reference nu values
for nu in range(len(matMassPropKeys)-1):
    ref_nu[nu] = ref_model.materials[nu+1].nu

# Load the mistuned model bdf
mistuned_model = BDF()
mistuned_model.read_bdf("inp/mistunedBeam.bdf", punch=True)
print(mistuned_model.get_bdf_stats())
   
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

for stiffness in range(len(matStiffPropKeys)-1):
    mistuned_stiffness[stiffness] = mistuned_model.materials[stiffness+1].e
    
# read mistuned mass values
for mass in range(len(matMassPropKeys)-1):
    mistuned_mass[mass] = mistuned_model.materials[mass+1].rho
    
# read reference nu values
for nu in range(len(matMassPropKeys)-1):
    mistuned_nu[nu] = mistuned_model.materials[nu+1].nu
    
    
# Prepare the arguments to be sent to SciPy using SQP

