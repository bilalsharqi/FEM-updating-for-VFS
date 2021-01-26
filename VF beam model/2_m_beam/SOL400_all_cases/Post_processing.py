

import numpy as np
import pandas as pd
import os
import shutil
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

plt.close("all")

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
file_path = ''
f06_file = 'beam_model_sol400.f06'
model_coords = 'sol400_coor.txt'
#Loads = np.asarray(pd.read_csv('Loads.txt', header=None))
n_modes=25
n_subcases=5

# BS - read the relevant files and store the numerical data
grids, n_grids, grid_coords = importGrids('', ['beam_model.bdf',model_coords], debug=True)
freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
# mode_shapes_NASTRAN = importEigenvectors(file_path, f06_file, n_modes, n_grids, grids, n_subcases,[],debug=True)
static_deform= importDisplacements(file_path, f06_file, n_subcases, grids, grids_order=[], debug=True)

print("\nNASTRAN data import completed")
    
fig=plt.figure(figsize=(15,8))
plt.plot(grid_coords[0,:] + static_deform[0][0][:], static_deform[0][2][:],'r',label='No load' )
plt.plot(grid_coords[0,:] + static_deform[1][0][:], static_deform[1][2][:],'b',label='OOP load')
plt.plot(grid_coords[0,:] + static_deform[2][0][:], static_deform[2][2][:],'k',label='IP load')
plt.plot(grid_coords[0,:] + static_deform[3][0][:], static_deform[3][2][:],'g',label='Rx load')
plt.plot(grid_coords[0,:] + static_deform[4][0][:], static_deform[4][2][:],'c',label='OOP + IP + Rx load')
plt.rcParams["lines.linewidth"] = 3
# plt.title('Static deformation with loading',fontsize=32 )
plt.xlabel("Beam length [m]",fontsize=26)
plt.ylabel("Vertical displacement [m]",fontsize=26)
plt.ax = plt.gca()
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=22)
plt.ylim(-0.1,1)

# fig.savefig("beam_static_2m_z.svg",bbox_inches='tight')
    
fig1=plt.figure(figsize=(15,8))
plt.plot(grid_coords[0,:] + static_deform[0][0][:], static_deform[0][1][:],'r',label='No load' )
plt.plot(grid_coords[0,:] + static_deform[1][0][:], static_deform[1][1][:],'b',label='OOP load')
plt.plot(grid_coords[0,:] + static_deform[2][0][:], static_deform[2][1][:],'k',label='IP load')
plt.plot(grid_coords[0,:] + static_deform[3][0][:], static_deform[3][1][:],'g',label='Rx load')
plt.plot(grid_coords[0,:] + static_deform[4][0][:], static_deform[4][1][:],'c',label='OOP + IP + Rx load')
plt.rcParams["lines.linewidth"] = 3
# plt.title('Static deformation with loading',fontsize=32 )
plt.xlabel("Beam length [m]",fontsize=26)
plt.ylabel("In-plane displacement [m]",fontsize=26)
plt.ax = plt.gca()
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=22)
plt.ylim(-0.1,.01)

# fig1.savefig("beam_static_2m_y.svg",bbox_inches='tight')

# BS - switch plot flag to yes if want to plot
# plot = 'yes'
# while plot=='yes':
    
    # Plotting
    # j=1 # subcase number
    # i=1 # Test mode number
    # print("MODE SHAPE PLOTTING")
    # print("OOP Bending mode shapes")
    
    # plt.figure()
    # plt.plot(grid_coords[0,:],grid_coords[2,:] + mode_shapes_NASTRAN[j][i][2][:],'k*',label='Z translation')
    # plt.title('Mode number ' + str(i) + ' at frequency ' + str(freq_NASTRAN[j][i]) + '' )
    # plt.xlabel('Length [m]')
    # plt.ylim(-1.5,1.5)
    # plt.ylabel('Displacement')
    # plt.legend()
    
    #print("IP and torsional components")    
    #plt.figure()
    #i=7 # Test mode number
    ##plt.plot(grid_coords[1,:],mode_shapes_NASTRAN[0][i][1][:],'r*',label='In-Plane component')
    #plt.plot(grid_coords[1,:],mode_shapes_NASTRAN[0][i][4][:],'b*',label='Rotation about y')
    #plt.title('Mode number ' + str(i) + ' at frequency ' + str(freq_NASTRAN[0][i]) + '' )
    #plt.legend()
    
    # Plot a basic scatter plot with static displacement 
    # plotted on top of the grid coordinates
    # print("STATIC DISPLACEMENT")
    # print("Translational data only")
    # fig = plt.figure(figsize=(16,9))
    # ax = fig.add_subplot(111, projection='3d')
    
    # ax.scatter(grid_coords[0,:] + static_deform[j][0][:], \
    #            grid_coords[1,:] + static_deform[j][1][:], \
    #            grid_coords[2,:] + static_deform[j][2][:])
    # plt.title('Static Displacement' )
    # #ax.set_xlim3d(-1,1.5)
    # #ax.set_ylim3d(-0.1,0.1)
    # ax.set_zlim3d(-1.5,1.5)
    # ax.set_xlabel('X [m]')
    # ax.set_ylabel('Y [m]')
    # ax.set_zlabel('Vertical displacement [m]')
    # #ax.invert_xaxis()
    # plt.show()
    
    # Plot a basic scatter plot with translational eigenvectors 
    # # plotted on top of the grid coordinates
    # print("MODE SHAPE PLOTTING")
    # print("Translational data only")
    # fig = plt.figure(figsize=(16,9))
    # ax = fig.add_subplot(111, projection='3d')
    # #i=6
    # ax.scatter(grid_coords[0,:] + mode_shapes_NASTRAN[j][i][0], \
    #            grid_coords[1,:] + mode_shapes_NASTRAN[j][i][1], \
    #            grid_coords[2,:] + mode_shapes_NASTRAN[j][i][2])
    # # the second scatter plot contains rotational data, which is not
    # # necessarily needed to visualize rotations
    # #ax.scatter(grid_coords[0,:]+mode_shapes_NASTRAN[0][i][3], \
    # #           grid_coords[1,:]+mode_shapes_NASTRAN[0][i][4], \
    # #           grid_coords[2,:]+mode_shapes_NASTRAN[0][i][5])
    
    # #ax.set_ylim3d(-0.1,0.1)
    # ax.set_zlim3d(-1.5,1.5)
    # ax.set_xlabel('X [m]')
    # ax.set_ylabel('Y [m]')
    # ax.set_zlabel('Vertical displacement [m]')
    # plt.title('Mode number ' + str(i) + ' at frequency ' + str(freq_NASTRAN[j][i]) + '' )
    # plt.show()
    
    # # Plot a basic scatter plot with translational eigenvectors 
    # # plotted on top of the grid coordinates which have the initial deformation
    # # coming from the static displacement accounted for
    # print("MODE SHAPE PLOTTING")
    # print("Deformed jig + Translational data only")
    # fig = plt.figure(figsize=(16,9))
    # ax = fig.add_subplot(111, projection='3d')
    # #i=6 # test mode shape number
    # # deformed mode shapes (initial grid coordinates + static displacement + modal displacements)
    # ax.scatter(grid_coords[0,:] + mode_shapes_NASTRAN[j][i][0] +  static_deform[j][0], \
    #            grid_coords[1,:] + mode_shapes_NASTRAN[j][i][1] +  static_deform[j][1], \
    #            grid_coords[2,:] + mode_shapes_NASTRAN[j][i][2] +  static_deform[j][2])
    # # static displacement (initial grid coordinates + static displacement)
    # ax.scatter(grid_coords[0,:] + static_deform[j][0][:], \
    #            grid_coords[1,:] + static_deform[j][1][:], \
    #            grid_coords[2,:] + static_deform[j][2][:])
    # # undeformed mode shapes (initial grid coordinates + modal displacements)
    # ax.scatter(grid_coords[0,:] + mode_shapes_NASTRAN[j][i][0], \
    #            grid_coords[1,:] + mode_shapes_NASTRAN[j][i][1], \
    #            grid_coords[2,:] + mode_shapes_NASTRAN[j][i][2])
    # #ax.set_ylim3d(-0.1,0.1)
    # ax.set_zlim3d(-1.5,1.5)
    # ax.set_xlabel('X [m]')
    # ax.set_ylabel('Y [m]')
    # ax.set_zlabel('Vertical displacement [m]')
    # #ax.invert_xaxis()
    # plt.title('Deformed Mode number ' + str(i) + ' at frequency ' + str(freq_NASTRAN[j][i]) + '' )
    # plt.show()
    
    # break

#from MTK.MTK import EigenPair, ModeSet, TrackModes
#
##================================================================================
## Mode tracking data creation: this is the formatting of data for MTK
##================================================================================
#
#print("\nStart of mode tracking data formatting")
#data_mode_sets =[[[0.,np.zeros([6*n_grids])] for i in range(n_modes)] for j in range(n_subcases)]
#for i in range(n_subcases):
#    for j in range(n_modes):
#        data_mode_sets[i][j][0] = freq_NASTRAN[i][j]
#        data_mode_sets[i][j][1] = mode_shapes_NASTRAN[i][j].flatten('F') # double check that it flattened correctly
#        
## function to plot data
#def PlotReal(var, data, line=True, sym=True):
#    """Plots the real mode progression over a variable.
#    """
#    N_sets = len(data)
#    N_modes = data[0].Size()
#    real = np.zeros([N_sets, N_modes])
#
#    for i in range(N_sets):
#        for j in range(N_modes):
#            real[i,j] = data[i][j]["value"].real
#
#    # plot the data
#    opt = ""
#    if sym:
#        opt += "o"
#    if line:
#        opt += "-"
#
#    plt.plot(var, real, opt)
#    
##print("\nRead and plot untracked data")
#
## create mode_sets to store each subcase as an entry in mode_sets
## NOTE: EigenPair(), ModeSet(), AddPair() are MTK functions
#mode_sets = []
#for modeset_list in data_mode_sets:
#    
#    modeset = ModeSet()
#    for mode in modeset_list:
#        pair = EigenPair()
#        pair["value"], pair["vector"]  = mode
#
#        modeset.AddPair(pair)
#
#    mode_sets.append(modeset)
#
## create seed or reference starting data    
#seed = ModeSet()
#
#for i in range(mode_sets[0].Size()):
##    print(i)
#    if (i >=0) and (i < 25):
##        print("adding mode " +str(i))
#        pair  = mode_sets[0][i]
#        seed.AddPair(pair)
#        
#var1=np.array(Loads)[0:,0] # store loads
## plot the untracked mode sets
#print("\nPlot untracked data")
#fig1=plt.figure(figsize=(15,8))
#PlotReal(var1, mode_sets[0:], sym=True, line=True)
#plt.rcParams["lines.linewidth"] = 3
#plt.title('Modes evolution with loading in g (untracked)',fontsize=32 )
#plt.xlabel("Loads [g]",fontsize=26)
#plt.ylabel("Frequency [Hz]",fontsize=26)
#plt.ax = plt.gca()
#plt.xticks(fontsize=25)
#plt.yticks(fontsize=25)
#plt.ax.spines["right"].set_visible(False)
#plt.ax.spines["top"].set_visible(False)
#
##fig1.savefig("untracked.svg",bbox_inches='tight')
#
##print("\ntrack modes")
## usage of TrackModes function in MTK
#tracked_mode_sets=TrackModes(seed,mode_sets[0:]) 
#var2=np.array(Loads)[0:,0]
#print("\nPlot tracked data")
#fig2=plt.figure(figsize=(15,8))
#plt.rcParams["lines.linewidth"] = 4
#PlotReal(var2, tracked_mode_sets, sym=True, line=True)
#plt.title('Modes evolution with loading: (tracked)',fontsize=32 )
#plt.xlabel("Loads [g]",fontsize=26)
#plt.xticks(fontsize=25)
#plt.yticks(fontsize=25)
#plt.ylabel("Frequency [Hz]",fontsize=26)
#plt.ax = plt.gca()
##plt.annotate('1 T', xy=(-0.05,46), xytext=(-0.05, 46),size=25)
##plt.annotate('2 IP', xy=(-0.05,53), xytext=(-0.05, 53),size=25)
##plt.annotate('2 T', xy=(-0.05,97), xytext=(-0.05, 97),size=25)
##plt.annotate('1 IP', xy=(-0.05, 18), xytext=(-0.05, 18),size=25)
##plt.annotate('Pitch:X', xy=(-0.05, 1), xytext=(-0.05, 1),size=25)
#plt.ax.spines["right"].set_visible(False)
#plt.ax.spines["top"].set_visible(False)
#
##fig2.savefig("tracked.svg",bbox_inches='tight')
#
### function to plot selected modes and save data in a .mat file if asked for
##def PlotReal2(var, data, line=True, sym=True):
##    """Plots the real mode progression over a variable. 
##    Modified to plot only certain modes input by the user in the 
##    test_modes array
##    """
##    N_sets = len(data)
##    N_modes = data[0].Size()
##    real = np.zeros([N_sets, N_modes])
##
##    for i in range(N_sets):
##        for j in range(N_modes):
##            real[i,j] = data[i][j]["value"].real
##
##    # plot the data
##    opt = ""
##    if sym:
##        opt += "o"
##    if line:
##        opt += "-"
##    test_modes=np.array([0,11,16,17,22])
##    real2=real[:,test_modes]
##    plt.plot(var, real2, opt)
##
##    dir=os.path.dirname(os.path.abspath("Post_processing.py"))
##    path = os.path.join(dir, "tracked_modes.mat")
##    database = {}
##    
##    # Write problem data
##    database["Loads"] = var2
##    database["tracked_modes"] = real2
##    
##    # Writing database, uncomment if want to store selected modes in a .mat file
###    print("...Exporting results in a .mat file")
###    if os.path.isfile(path):
###        os.remove(path)
###    sio.savemat(path,database,appendmat=False)
##    
###print("\nselected tracked modes")
##tracked_mode_sets=TrackModes(seed,mode_sets[0:])
##var2=np.array(Loads)[0:,0]
##    
##print("\nPlot tracked data with modes of interest")
##fig4=plt.figure(figsize=(15,8))
##plt.rcParams["lines.linewidth"] = 4
##PlotReal2(var2, tracked_mode_sets, sym=True, line=True)
##plt.title('Selected modes evolution with loading: (tracked)',fontsize=32 )
##plt.xlabel("Loads [g]",fontsize=26)
##plt.xticks(fontsize=25)
##plt.yticks(fontsize=25)
##plt.ylabel("Frequency [Hz]",fontsize=26)
##plt.ax = plt.gca()
###plt.annotate('Pitch:X', xy=(-0.05, 1), xytext=(-0.05, 1),size=25)
###plt.annotate('1 IP', xy=(-0.05, 18),   xytext=(-0.05, 18),size=25)
###plt.annotate('1 T', xy=(-0.05,45),     xytext=(-0.05, 45),size=25)
###plt.annotate('2 IP', xy=(-0.05,53),    xytext=(-0.05, 53),size=25)
###plt.annotate('2 T', xy=(-0.05,97),     xytext=(-0.05, 97),size=25)
##
##plt.ax.spines["right"].set_visible(False)
##plt.ax.spines["top"].set_visible(False)
##
###fig4.savefig("tracked_test_out_tuned_selected.svg",bbox_inches='tight')