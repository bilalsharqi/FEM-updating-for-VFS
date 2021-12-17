# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 19:11:03 2021

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
from pyMSCNastranUtilities import *
plt.close('all')
# import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
temp_path = 'Results_11_12_21_NM_GF_UNCON/'
# temp_path = 'Results_11_12_21_SLSQP_BG_CON_SOL103/'
# temp_path = 'Results_11_15_21_SLSQP_GB_CON/'
# temp_path = 'Results_11_28_21_SLSQP_GB_CON_SOL400_w_static/'
n_modes=15
n_subcases=4

# Read reference results
# Note: There are 4 loading subcases and 15 eigenvalues computed for each
# deformed modal analysis
file_path = temp_path + 'out'
f06_file = 'sol400.f06'
model_coords = 'sol400_coor.txt'

# Components

# read the reference result files and store the numerical data
ref_grids, ref_n_grids, ref_grid_coords = importGrids(file_path, ['refBeam.bdf',model_coords], debug=True)
ref_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
ref_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, ref_n_grids, ref_grids, n_subcases,[],debug=True)
ref_static_deform= importDisplacements(file_path, f06_file, n_subcases, ref_grids, grids_order=[], debug=True)
M_ref, x_G_ref, J_G_ref = importRigidBodyMassData(file_path, f06_file,debug=True)
print("\nReference NASTRAN data import completed")

# Read initial mistuned results
# Note: There are 4 loading subcases and 15 eigenvalues computed for each
# deformed modal analysis
file_path = temp_path + 'inp'
f06_file = 'mistuned_sol400.f06'
model_coords = 'mistuned_sol400_coor.txt'

# read the reference result files and store the numerical data
initial_mistuned_grids, initial_mistuned_n_grids, initial_mistuned_grid_coords = importGrids(file_path, ['mistunedBeam.bdf',model_coords], debug=True)
initial_mistuned_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
initial_mistuned_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, initial_mistuned_n_grids, initial_mistuned_grids, n_subcases,[],debug=True)
initial_mistuned_static_deform= importDisplacements(file_path, f06_file, n_subcases, initial_mistuned_grids, grids_order=[], debug=True)
M_mistuned, x_G_mistuned, J_G_mistuned = importRigidBodyMassData(file_path, f06_file,debug=True)
print("\nInitial mistuned NASTRAN data import completed")


# Read final/converged mistuned results
# Note: There are 4 loading subcases and 15 eigenvalues computed for each
# deformed modal analysis
file_path = temp_path + 'out'
f06_file = 'mistuned_sol400.f06'
model_coords = 'mistuned_sol400_coor.txt'

# read the reference result files and store the numerical data
mistuned_grids, mistuned_n_grids, mistuned_grid_coords = importGrids(file_path, ['mistunedBeam1.bdf',model_coords], debug=True)
mistuned_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
mistuned_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, mistuned_n_grids, mistuned_grids, n_subcases,[],debug=True)
mistuned_static_deform= importDisplacements(file_path, f06_file, n_subcases, mistuned_grids, grids_order=[], debug=True)
M_mistuned, x_G_mistuned, J_G_mistuned = importRigidBodyMassData(file_path, f06_file,debug=True)
print("\nMistuned NASTRAN data import completed")


subcase = 0

# Static deformation
fig=plt.figure(figsize=(17,10))
plt.plot(ref_grid_coords[0,:] + ref_static_deform[0][0][:]-1, ref_static_deform[0][2][:]/1.5*100,'kx',label='Reference FEM: -1 g load' ,markersize=8)
plt.plot(ref_grid_coords[0,:] + ref_static_deform[1][0][:]-1, ref_static_deform[1][2][:]/1.5*100,'k^',label='Reference FEM: 10 N distributed load',markersize=8)
plt.plot(ref_grid_coords[0,:] + ref_static_deform[2][0][:]-1, ref_static_deform[2][2][:]/1.5*100,'ko',label='Reference FEM: 25 N distributed load',markersize=8)
plt.plot(ref_grid_coords[0,:]-1, ref_grid_coords[2,:]-1,'k+',label='Reference FEM: No load',markersize=8)
plt.plot(initial_mistuned_grid_coords[0,:] + initial_mistuned_static_deform[0][0][:]-1, initial_mistuned_static_deform[0][2][:]/1.5*100,label='Initial FEM: -1 g load' ,markersize=6)
plt.plot(initial_mistuned_grid_coords[0,:] + initial_mistuned_static_deform[1][0][:]-1, initial_mistuned_static_deform[1][2][:]/1.5*100,label='Initial FEM: 10 N distributed load',markersize=6)
plt.plot(initial_mistuned_grid_coords[0,:] + initial_mistuned_static_deform[2][0][:]-1, initial_mistuned_static_deform[2][2][:]/1.5*100,label='Initial FEM: 25 N distributed load',markersize=6)
plt.plot(initial_mistuned_grid_coords[0,:]-1, initial_mistuned_grid_coords[2,:]-1,label='Initial FEM: No load',markersize=6)
plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[0][0][:]-1, mistuned_static_deform[0][2][:]/1.5*100,label='Final FEM: -1 g load' ,markersize=4)
plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[1][0][:]-1, mistuned_static_deform[1][2][:]/1.5*100,label='Final FEM: 10 N distributed load',markersize=4)
plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[2][0][:]-1, mistuned_static_deform[2][2][:]/1.5*100,label='Final FEM: 25 N distributed load',markersize=4)
plt.plot(mistuned_grid_coords[0,:]-1, mistuned_grid_coords[2,:]-1,label='Final FEM: No load',markersize=4)
plt.rcParams["lines.linewidth"] = 3
# plt.title('Static deformation with loading',fontsize=32 )
plt.xlabel("Beam length [m]",fontsize=26)
plt.ylabel("Vertical displacement [% span]",fontsize=26)
plt.ax = plt.gca()
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=16)
plt.ax.spines['top'].set_visible(False)
plt.ax.spines['right'].set_visible(False)
# plt.ax.spines['bottom'].set_visible(False)
# plt.ax.spines['left'].set_visible(False)
# plt.ylim(-1,1.5)
fig.savefig("static_displacement/beam_static_displacement.svg",bbox_inches='tight')


# print("Deformed OOP Bending mode shapes")
# # select component of mode shape (3 translations and 3 rotations)
# phi_component = 2
# # for i in range(0,n_modes):
# i=3 # for quickly plotting any one mode
# fig=plt.figure(figsize=(15,8))
# plt.plot(ref_grid_coords[0,:]-1,ref_mode_shapes[subcase][i][phi_component]+ref_static_deform[subcase][phi_component],'k^',label='Reference FEM')
# plt.plot(initial_mistuned_grid_coords[0,:]-1,initial_mistuned_mode_shapes[subcase][i][phi_component]+initial_mistuned_static_deform[subcase][phi_component],'b',label='Initial mistuned FEM')
# plt.plot(mistuned_grid_coords[0,:]-1,mistuned_mode_shapes[subcase][i][phi_component]+mistuned_static_deform[subcase][phi_component],'r',label='Converged FEM')
# # plt.title('Mode number ' + str(i) + ' at frequency ' + str(ref_freq_NASTRAN[subcase][i]) + '' ,fontsize=26)
# plt.xlabel('Beam Length [m]',fontsize=26)
# plt.rcParams["lines.linewidth"] = 3
# plt.ylabel('Amplitude',fontsize=26)
# plt.ax = plt.gca()
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# plt.legend(fontsize=22)
# plt.ax.spines['top'].set_visible(False)
# plt.ax.spines['right'].set_visible(False)
# # plt.ax.spines['bottom'].set_visible(False)
# # plt.ax.spines['left'].set_visible(False)
# fig.savefig('mode_shapes/Mode shape ' + str(i+1) + '.svg',bbox_inches='tight')

# Plot frequencies
# n_freq = 5
# for i in range(0, n_subcases):
#     fig=plt.figure(figsize=(16,9))
#     plt.plot(np.linspace(1,n_freq,n_freq),ref_freq_NASTRAN[i][0:5],'k^',label='Reference FEM',markersize=8)
#     plt.plot(np.linspace(1,n_freq,n_freq),mistuned_freq_NASTRAN[i][0:5],'ro',label='Converged FEM',markersize=8)
#     plt.plot(np.linspace(1,n_freq,n_freq),initial_mistuned_freq_NASTRAN[i][0:5],'b+',label='Initial FEM',markersize=12)
#     plt.xlabel('Mode number',fontsize=26)
#     plt.rcParams["lines.linewidth"] = 3
#     plt.ylabel('Frequency [Hz]',fontsize=26)
#     plt.ax = plt.gca()
#     plt.xticks(fontsize=25)
#     plt.ylim(0,40)
#     plt.xticks(np.arange(1, n_freq+1, step=1.0))
#     plt.yticks(fontsize=25)
#     plt.legend(fontsize=22)
#     plt.ax.spines['top'].set_visible(False)
#     plt.ax.spines['right'].set_visible(False)
#     fig.savefig('mode_shapes/Frequency case ' + str(i+1) + '.svg',bbox_inches='tight')
    
# # MAC module
# temp_MAC_size = (n_subcases,n_subcases)

# MAC_matrix_ref_initial = [np.zeros(temp_MAC_size) for j in range(n_subcases)]
# MAC_matrix_ref_final = [np.zeros(temp_MAC_size) for j in range(n_subcases)]

# data_MAC_reference = [[np.zeros([6*mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]
# data_MAC_initial  = [[np.zeros([6*initial_mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]
# data_MAC_mistuned  = [[np.zeros([6*mistuned_n_grids]) for i in range(n_modes)] for j in range(n_subcases)]

# for i in range(n_subcases):
#     for j in range(n_modes):
#         data_MAC_reference[i][j] = ref_mode_shapes[i][j].flatten('C')
#         data_MAC_initial[i][j] = initial_mistuned_mode_shapes[i][j].flatten('C')
#         data_MAC_mistuned[i][j] = mistuned_mode_shapes[i][j].flatten('C')

# # test_MAC2 = ComputeMAC(data_MAC_reference[0][0],data_MAC_mistuned[0][0])
# for i in range(n_subcases):
#         MAC_matrix_ref_initial[i] = mac(np.array(data_MAC_reference[i][:]), np.array(data_MAC_initial[i][:]))
#         MAC_matrix_ref_final[i] = mac(np.array(data_MAC_reference[i][:]), np.array(data_MAC_mistuned[i][:]))
        
# # Plot initial MAC        
# fig = plt.figure(figsize=(16, 9))
# ax = fig.add_subplot(111)
# ax.set_aspect('equal')
# im=plt.imshow(MAC_matrix_ref_initial[subcase])
# cb=plt.colorbar(im,orientation='vertical')
# cb.ax.tick_params(labelsize=20)
# plt.clim(0, 1.0) 
# im.figure.axes[0].tick_params(axis="both", labelsize=20)
# im.figure.axes[1].tick_params(axis="x", labelsize=20)
# plt.show()
# fig.savefig("MAC/initial_vs_ref_MAC_case_0.svg",bbox_inches='tight')

# # Plot final MAC        
# fig = plt.figure(figsize=(16, 9))
# ax = fig.add_subplot(111)
# ax.set_aspect('equal')
# im=plt.imshow(MAC_matrix_ref_final[subcase])
# cb=plt.colorbar(im,orientation='vertical')
# cb.ax.tick_params(labelsize=20)
# plt.clim(0, 1.0) 
# im.figure.axes[0].tick_params(axis="both", labelsize=20)
# im.figure.axes[1].tick_params(axis="x", labelsize=20)
# plt.show()
# fig.savefig("MAC/final_vs_ref_MAC_case_0.svg",bbox_inches='tight')

