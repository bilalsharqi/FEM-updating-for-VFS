

import numpy as np
import pandas as pd
import os
import shutil
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

plt.close("all")

sys.path.append('../Individual Load Cases/1 g/')
sys.path.append('../Individual Load Cases/1 g/')

from pyMSCNastranUtilities import *

# BS - if a single load case is desired, uncomment lines below
#file_path='../Individual Load Cases/0 g/'
#file_path = '../Individual Load Cases/1 g/'
#f06_file = 'beam_model_sol400.f06'
#model_coords = 'sol400_coor.txt'
#n_modes=25
#n_subcases=1

# BS - if multiple load cases are desired, import these files 
# NOTE: I have not finished collating this dataset for all loads yet
file_path = '../Multiple Load Cases/'
f06_file = 'beam_model_sol400_case_1to48.f06'
model_coords = 'sol400_coor.txt'
Loads = np.asarray(pd.read_csv('Loads.txt'))
n_modes=25
n_subcases=48

# BS - read the relevant files and store the numerical data
grids, n_grids, grid_coords = importGrids('', ['beam_model.bdf',model_coords], debug=True)
freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
mode_shapes_NASTRAN = importEigenvectors(file_path, f06_file, n_modes, n_grids, grids, n_subcases,[],debug=True)
static_deform= importDisplacements(file_path, f06_file, n_subcases, grids, grids_order=[], debug=True)

print("\nNASTRAN data import completed")
    
# BS - switch plot flag to yes if want to plot
plot = 'no'
while plot=='yes':
    
    # Plotting
    j=47 # subcase number
    i=7 # Test mode number
    print("MODE SHAPE PLOTTING")
    print("OOP Bending mode shapes")
    
    plt.figure()
    plt.plot(grid_coords[0,:],grid_coords[2,:] + mode_shapes_NASTRAN[j][i][2][:],'k*',label='Z translation')
    plt.title('Mode number ' + str(i) + ' at frequency ' + str(freq_NASTRAN[j][i]) + '' )
    plt.xlabel('Length [m]')
    plt.ylim(-1.5,1.5)
    plt.ylabel('Displacement')
    plt.legend()
    
    #print("IP and torsional components")    
    #plt.figure()
    #i=7 # Test mode number
    ##plt.plot(grid_coords[1,:],mode_shapes_NASTRAN[0][i][1][:],'r*',label='In-Plane component')
    #plt.plot(grid_coords[1,:],mode_shapes_NASTRAN[0][i][4][:],'b*',label='Rotation about y')
    #plt.title('Mode number ' + str(i) + ' at frequency ' + str(freq_NASTRAN[0][i]) + '' )
    #plt.legend()
    
    # Plot a basic scatter plot with static displacement 
    # plotted on top of the grid coordinates
    print("STATIC DISPLACEMENT")
    print("Translational data only")
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.scatter(grid_coords[0,:] + static_deform[j][0][:], \
               grid_coords[1,:] + static_deform[j][1][:], \
               grid_coords[2,:] + static_deform[j][2][:])
    plt.title('Static Displacement' )
    #ax.set_xlim3d(-1,1.5)
    #ax.set_ylim3d(-0.1,0.1)
    ax.set_zlim3d(-1.5,1.5)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Vertical displacement [m]')
    #ax.invert_xaxis()
    plt.show()
    
    # Plot a basic scatter plot with translational eigenvectors 
    # plotted on top of the grid coordinates
    print("MODE SHAPE PLOTTING")
    print("Translational data only")
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111, projection='3d')
    #i=6
    ax.scatter(grid_coords[0,:] + mode_shapes_NASTRAN[j][i][0], \
               grid_coords[1,:] + mode_shapes_NASTRAN[j][i][1], \
               grid_coords[2,:] + mode_shapes_NASTRAN[j][i][2])
    # the second scatter plot contains rotational data, which is not
    # necessarily needed to visualize rotations
    #ax.scatter(grid_coords[0,:]+mode_shapes_NASTRAN[0][i][3], \
    #           grid_coords[1,:]+mode_shapes_NASTRAN[0][i][4], \
    #           grid_coords[2,:]+mode_shapes_NASTRAN[0][i][5])
    
    #ax.set_ylim3d(-0.1,0.1)
    ax.set_zlim3d(-1.5,1.5)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Vertical displacement [m]')
    plt.title('Mode number ' + str(i) + ' at frequency ' + str(freq_NASTRAN[j][i]) + '' )
    plt.show()
    
    # Plot a basic scatter plot with translational eigenvectors 
    # plotted on top of the grid coordinates which have the initial deformation
    # coming from the static displacement accounted for
    print("MODE SHAPE PLOTTING")
    print("Deformed jig + Translational data only")
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111, projection='3d')
    #i=6 # test mode shape number
    # deformed mode shapes (initial grid coordinates + static displacement + modal displacements)
    ax.scatter(grid_coords[0,:] + mode_shapes_NASTRAN[j][i][0] +  static_deform[j][0], \
               grid_coords[1,:] + mode_shapes_NASTRAN[j][i][1] +  static_deform[j][1], \
               grid_coords[2,:] + mode_shapes_NASTRAN[j][i][2] +  static_deform[j][2])
    # static displacement (initial grid coordinates + static displacement)
    ax.scatter(grid_coords[0,:] + static_deform[j][0][:], \
               grid_coords[1,:] + static_deform[j][1][:], \
               grid_coords[2,:] + static_deform[j][2][:])
    # undeformed mode shapes (initial grid coordinates + modal displacements)
    ax.scatter(grid_coords[0,:] + mode_shapes_NASTRAN[j][i][0], \
               grid_coords[1,:] + mode_shapes_NASTRAN[j][i][1], \
               grid_coords[2,:] + mode_shapes_NASTRAN[j][i][2])
    #ax.set_ylim3d(-0.1,0.1)
    ax.set_zlim3d(-1.5,1.5)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Vertical displacement [m]')
    #ax.invert_xaxis()
    plt.title('Deformed Mode number ' + str(i) + ' at frequency ' + str(freq_NASTRAN[j][i]) + '' )
    plt.show()
    
    break
