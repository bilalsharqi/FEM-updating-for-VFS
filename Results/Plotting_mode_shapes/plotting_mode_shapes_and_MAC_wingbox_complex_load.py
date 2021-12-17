# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 19:11:03 2021

@author: bilal
"""

import numpy as np
import os
import scipy.optimize as sio
import scipy.io as s_io
import shutil
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib import cm
from pyMSCNastranUtilities import *
plt.close('all')

temp_path = 'Renato_wingbox/'
n_modes=15
n_subcases=1


# Read final/converged mistuned results
# Note: There are 4 loading subcases and 15 eigenvalues computed for each
# deformed modal analysis
file_path = temp_path + 'reference/Complex_loads/Renumbered_grids'
f06_file = 'sol400_complex_loads.f06'
model_coords = 'mistuned_sol400_coor_cplx_load.txt'

# read the reference result files and store the numerical data
mistuned_grids, mistuned_n_grids, mistuned_grid_coords = importGrids(file_path, ['wingbox.bdf',model_coords], debug=True)
mistuned_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
mistuned_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, mistuned_n_grids, mistuned_grids, n_subcases,[],debug=True)
mistuned_static_deform= importDisplacements(file_path, f06_file, n_subcases, mistuned_grids, grids_order=[], debug=True)
M_mistuned, x_G_mistuned, J_G_mistuned = importRigidBodyMassData(file_path, f06_file,debug=True)
print("\nMistuned NASTRAN data import completed")

# # ================================================================================
# # Save data
# # ================================================================================
print("...Exporting results in a .mat file")
dir=os.path.dirname(os.path.abspath("plotting_mode_shapes_and_MAC_wingbox.py"))
path = os.path.join(dir, "visualize_complex_load_reference.mat")
database = {}

# Write problem data
database["n_subcases"] = n_subcases
database['n_modes'] = n_modes
database["n_grids"] = mistuned_n_grids
database["ref_grid_coords"] = mistuned_grids

# Write mode shapes
database['mode_shapes'] = mistuned_mode_shapes
# Write frequencies 
database["frequencies"] = mistuned_freq_NASTRAN
# write static displacement
database['static_displacement'] = mistuned_static_deform

# Writing database
if os.path.isfile(path):
    os.remove(path)
s_io.savemat(path,database,appendmat=False)


# subcase = 0

# # Plot a basic scatter plot with static displacement 
# # plotted on top of the grid coordinates
# print("STATIC DISPLACEMENT")
# print("Translational data only")
# fig = plt.figure(figsize=(16,9))
# ax = fig.add_subplot(111, projection='3d')
# ax.set_box_aspect((np.ptp(mistuned_grid_coords[0,:] + mistuned_static_deform[subcase][0]), 3*np.ptp(mistuned_grid_coords[1,:] + mistuned_static_deform[subcase][1]), np.ptp(mistuned_grid_coords[2,:] + mistuned_static_deform[subcase][2])))

# # ax.scatter3D(ref_grid_coords[0,:] + ref_static_deform[subcase][0], \
# #            ref_grid_coords[1,:] + ref_static_deform[subcase][1], \
# #            ref_grid_coords[2,:] + ref_static_deform[subcase][2],s=50,marker='.',c='k',label='Reference FEM')
# # ax.scatter3D(initial_mistuned_grid_coords[0,:] + initial_mistuned_static_deform[subcase][0], \
# #            initial_mistuned_grid_coords[1,:] + initial_mistuned_static_deform[subcase][1], \
# #            initial_mistuned_grid_coords[2,:] + initial_mistuned_static_deform[subcase][2],s=1,marker='.',c='b',label='Initial FEM')
# ax.scatter3D(mistuned_grid_coords[0,:] + mistuned_static_deform[subcase][0], \
#            mistuned_grid_coords[1,:] + mistuned_static_deform[subcase][1], \
#            mistuned_grid_coords[2,:] + mistuned_static_deform[subcase][2],s=1,marker='.',c='r',label='Converged FEM')
# plt.title('Static Displacement' )
# plt.ax = plt.gca()
# # plt.xticks(fontsize=25)
# # plt.yticks(fontsize=25)
# # plt.zticks(fontsize=25)
# # plt.legend(fontsize=14)
# lgnd = plt.legend(loc="best", scatterpoints=1, fontsize=14)
# lgnd.legendHandles[0]._sizes = [60]
# # lgnd.legendHandles[1]._sizes = [60]
# # lgnd.legendHandles[2]._sizes = [60]
# #ax.set_xlim3d(-1,1.5)
# # ax.set_ylim3d(-3,3)
# # ax.set_zlim3d(-1.5,1.5)
# ax.set_xlabel('Span [m]')
# ax.set_ylabel('Chord [m]')
# ax.set_zlabel('Vertical displacement [m]')
# #ax.invert_xaxis()
# plt.show()


# # Plot a basic scatter plot with translational eigenvectors 
# # plotted on top of the grid coordinates which have the initial deformation
# # coming from the static displacement accounted for
# print("MODE SHAPE PLOTTING")
# print("Deformed jig + Translational data only")
# fig = plt.figure(figsize=(16,9))
# ax = fig.add_subplot(111, projection='3d')
# ax.set_box_aspect((np.ptp(ref_grid_coords[0,:] + ref_static_deform[subcase][0]), np.ptp(ref_grid_coords[1,:] + ref_static_deform[subcase][1]), np.ptp(ref_grid_coords[2,:] + ref_static_deform[subcase][2])))
# i=0 # test mode shape number
# # deformed mode shapes (initial grid coordinates + static displacement + modal displacements)
# # ax.scatter(ref_grid_coords[0,:] + ref_mode_shapes[0][i][0] +  ref_static_deform[0][0], \
# #            ref_grid_coords[1,:] + ref_mode_shapes[0][i][1] +  ref_static_deform[0][1], \
# #            ref_grid_coords[2,:] + ref_mode_shapes[0][i][2] +  ref_static_deform[0][2],s=1,marker='.')
# # # static displacement (initial grid coordinates + static displacement)
# # ax.scatter(ref_grid_coords[0,:] + ref_static_deform[0][0][:], \
# #            ref_grid_coords[1,:] + ref_static_deform[0][1][:], \
# #            ref_grid_coords[2,:] + ref_static_deform[0][2][:],s=1,marker='.')
# # # undeformed mode shapes (initial grid coordinates + modal displacements)
# ax.scatter(ref_grid_coords[0,:] + ref_mode_shapes[0][i][0], \
#             ref_grid_coords[1,:] + ref_mode_shapes[0][i][1], \
#             ref_grid_coords[2,:] + ref_mode_shapes[0][i][2],s=1,marker='.')
# # ax.set_ylim3d(-3,3)
# # ax.set_zlim3d(-1.5,1.5)
# ax.set_xlabel('Span [m]')
# ax.set_ylabel('Chord [m]')
# ax.set_zlabel('Vertical displacement [m]')
# #ax.invert_xaxis()
# plt.title('Deformed Mode number ' + str(i+1) + ' at frequency ' + str(ref_freq_NASTRAN[0][i]) + '' )
# plt.show()
# fig.savefig("Renato_wingbox/static_displacement/static_displacement8.svg",bbox_inches='tight')

# # Static deformation
# fig=plt.figure(figsize=(15,8))
# plt.plot(ref_grid_coords[0,:] + ref_static_deform[0][0][:]-1, ref_static_deform[0][2][:],'kx',label='Gravity load ref' ,markersize=8)
# plt.plot(ref_grid_coords[0,:] + ref_static_deform[1][0][:]-1, ref_static_deform[1][2][:],'k^',label='10 N distributed load ref',markersize=8)
# plt.plot(ref_grid_coords[0,:] + ref_static_deform[2][0][:]-1, ref_static_deform[2][2][:],'ko',label='25 N distributed load ref',markersize=8)
# plt.plot(ref_grid_coords[0,:]-1, ref_grid_coords[2,:]-1,'k+',label='No load ref',markersize=8)
# plt.plot(initial_mistuned_grid_coords[0,:] + initial_mistuned_static_deform[0][0][:]-1, initial_mistuned_static_deform[0][2][:],label='Gravity load ini' ,markersize=6)
# plt.plot(initial_mistuned_grid_coords[0,:] + initial_mistuned_static_deform[1][0][:]-1, initial_mistuned_static_deform[1][2][:],label='10 N distributed load ini',markersize=6)
# plt.plot(initial_mistuned_grid_coords[0,:] + initial_mistuned_static_deform[2][0][:]-1, initial_mistuned_static_deform[2][2][:],label='25 N distributed load ini',markersize=6)
# plt.plot(initial_mistuned_grid_coords[0,:]-1, initial_mistuned_grid_coords[2,:]-1,label='No load ini',markersize=6)
# plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[0][0][:]-1, mistuned_static_deform[0][2][:],label='Gravity load fin' ,markersize=4)
# plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[1][0][:]-1, mistuned_static_deform[1][2][:],label='10 N distributed load fin',markersize=4)
# plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[2][0][:]-1, mistuned_static_deform[2][2][:],label='25 N distributed load fin',markersize=4)
# plt.plot(mistuned_grid_coords[0,:]-1, mistuned_grid_coords[2,:]-1,label='No load fin',markersize=4)
# # plt.plot(grid_coords[0,:] + static_deform[3][0][:], static_deform[3][2][:],'g',label='No load')
# # plt.plot(grid_coords[0,:] + static_deform[4][0][:], static_deform[4][2][:],'c',label='OOP + IP + Rx load')
# plt.rcParams["lines.linewidth"] = 3
# # plt.title('Static deformation with loading',fontsize=32 )
# plt.xlabel("Beam length [m]",fontsize=26)
# plt.ylabel("Vertical displacement [m]",fontsize=26)
# plt.ax = plt.gca()
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# plt.legend(fontsize=14)
# # plt.ylim(-1,1.5)
# fig.savefig("Renato_wingbox/static_displacement/static_displacement.svg",bbox_inches='tight')


# print("Deformed OOP Bending mode shapes")
# # select component of mode shape (3 translations and 3 rotations)
# phi_component = 2
# for i in range(0,n_modes):
#     # i=2 for quickly plotting any one mode
#     fig=plt.figure(figsize=(15,8))
#     plt.plot(ref_grid_coords[0,:]-1,ref_mode_shapes[subcase][i][phi_component]+ref_static_deform[subcase][phi_component],'k^',label='Reference FEM')
#     plt.plot(initial_mistuned_grid_coords[0,:]-1,initial_mistuned_mode_shapes[subcase][i][phi_component]+initial_mistuned_static_deform[subcase][phi_component],'b',label='Initial mistuned FEM')
#     plt.plot(mistuned_grid_coords[0,:]-1,mistuned_mode_shapes[subcase][i][phi_component]+mistuned_static_deform[subcase][phi_component],'r',label='Converged FEM')
#     # plt.title('Mode number ' + str(i) + ' at frequency ' + str(ref_freq_NASTRAN[subcase][i]) + '' ,fontsize=26)
#     plt.xlabel('Beam Length [m]',fontsize=26)
#     plt.rcParams["lines.linewidth"] = 3
#     plt.ylabel('Amplitude',fontsize=26)
#     plt.ax = plt.gca()
#     plt.xticks(fontsize=25)
#     plt.yticks(fontsize=25)
#     plt.legend(fontsize=22)
#     fig.savefig('Renato_wingbox/mode_shapes/Mode shape ' + str(i+1) + '.svg',bbox_inches='tight')

# # Plot frequencies
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
#     fig.savefig('Renato_wingbox/mode_shapes/Frequency case ' + str(i+1) + '.svg',bbox_inches='tight')
    
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
# fig.savefig("Renato_wingbox/MAC/initial_vs_ref_MAC.svg",bbox_inches='tight')

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
# fig.savefig("Renato_wingbox/MAC/final_vs_ref_MAC.svg",bbox_inches='tight')



data = np.array([[1,	0.68,	0.68,	0.83,	0.68,	0.23,],
        [2,	2.87,	2.87,	0.03,	2.88,	0.10,],
        [3,	4.25,	4.27,	0.42,	4.24,	0.21,],
        [4,	11.71,	11.73,	0.24,	11.69,	0.16,],
        [5,	17.20,	17.20,	0.04,	17.20,	0.04,],
        [6,	18.04,	18.03,	0.05,	18.05,	0.02,],
        [7,	22.18,	22.22,	0.17,	22.16,	0.11,],
        [8,	34.81,	34.91,	0.26,	34.80,	0.04,],
        [9,	47.90,	47.87,	0.05,	47.93,	0.07,],
        [10,	48.23,	48.58,	0.72,	48.25,	0.05,],
        [11,	52.23,	52.24,	0.01,	52.18,	0.09,],
        [12,	54.89,	58.01,	5.69,	55.07,	0.33,],
        [13,	60.57,	62.54,	3.26,	60.64,	0.12,],
        [14,	63.92,	64.00,	0.12,	63.88,	0.06,],
        [15,	67.16,	69.97,	4.19,	67.25,	0.14]])

# data = np.array([[1	   , 0.51	,0.52,	0.63,	0.51	,0.16],
#                  [2	   , 2.71	,2.71,	0.11,	2.71	,0.19],
#                  [3	   , 3.90	,3.90,	0.14,	3.89	,0.09],
#                  [4	    ,10.55	,10.60	,0.47	,10.55	,0.03],
#                  [5	    ,16.78	,16.78	,0.01	,16.78	,0.01],
#                  [6	    ,17.80	,17.81	,0.03	,17.81	,0.02],
#                  [7	    ,19.97	,20.00	,0.19	,19.97	,0.00],
#                  [8	    ,31.77	,31.87	,0.32	,31.79	,0.05],
#                  [9	    ,43.50	,43.80	,0.68	,43.56	,0.15],
#                  [10	,47.07	,47.05	,0.04	,47.12	,0.09],
#                  [11	,51.43	,51.73	,0.57	,51.65	,0.42],
#                  [12	,51.69	,53.45	,3.40	,51.73	,0.07],
#                  [13	,56.26	,56.27	,0.01	,56.38	,0.21],
#                  [14	,59.09	,60.08	,1.67	,59.19	,0.18],
#                  [15	,60.81	,62.62	,2.98	,60.92	,0.18]])


# Plot errors for complex load case linear vs nonlinear optimization
n_freq = 15
for i in range(0, n_subcases):
    fig=plt.figure(figsize=(16,9))
    plt.plot(np.linspace(1,n_freq,n_freq),data[:,3],'k^',label='Conventional optimization vs. reference FEM',markersize=8)
    plt.plot(np.linspace(1,n_freq,n_freq),data[:,5],'ro',label='Modified optimization vs. reference FEM',markersize=8)
    plt.xlabel('Mode number',fontsize=26)
    plt.rcParams["lines.linewidth"] = 3
    plt.ylabel('Error [%]',fontsize=26)
    plt.ax = plt.gca()
    plt.xticks(fontsize=25)
    plt.ylim(0,10)
    plt.xticks(np.arange(1, n_freq+1, step=1.0))
    plt.yticks(fontsize=25)
    plt.legend(fontsize=22)
    plt.ax.spines['top'].set_visible(False)
    plt.ax.spines['right'].set_visible(False)
    fig.savefig('mode_shapes/linear_v_nonlinear.svg',bbox_inches='tight')
