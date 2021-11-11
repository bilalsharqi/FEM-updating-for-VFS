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
# runNastran(inp_dir, run_dir, out_dir, "sol400.dat", debug=True)

# Read reference results
file_path = "out"
f06_file = 'sol400.f06'
model_coords = 'sol400_coor.txt'
n_modes=15
n_subcases=4

# read the reference result files and store the numerical data
ref_grids, ref_n_grids, ref_grid_coords = importGrids(file_path, ['wingbox.bdf',model_coords], debug=True)
ref_freq_NASTRAN = importFrequencies(file_path, f06_file, n_modes, n_subcases, debug=True)
ref_mode_shapes = importEigenvectors(file_path, f06_file, n_modes, ref_n_grids, ref_grids, n_subcases,[],debug=True)
ref_static_deform= importDisplacements(file_path, f06_file, n_subcases, ref_grids, grids_order=[], debug=True)

print("\nNASTRAN data import completed")
    
fig=plt.figure(figsize=(15,8))
plt.plot(ref_grid_coords[0,:] + ref_static_deform[0][0][:]-1, ref_static_deform[0][2][:],'r',label='Gravity load' )
plt.plot(ref_grid_coords[0,:] + ref_static_deform[1][0][:]-1, ref_static_deform[1][2][:],'b',label='10 N distributed load')
plt.plot(ref_grid_coords[0,:] + ref_static_deform[2][0][:]-1, ref_static_deform[2][2][:],'g',label='25 N distributed load')
plt.plot(ref_grid_coords[0,:]-1, ref_grid_coords[2,:]-1, 'k',label='No load')
plt.rcParams["lines.linewidth"] = 3
plt.title('Reference Static deformation with loading',fontsize=32 )
plt.xlabel("Beam length [m]",fontsize=26)
plt.ylabel("Vertical displacement [m]",fontsize=26)
plt.ax = plt.gca()
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=22)
# plt.ylim(-1,1.5)

fig.savefig("out/beam_static_1_5m_z.svg",bbox_inches='tight')

# Call Nastran for generating mistuned results
# runNastran(inp_dir, run_dir, out_dir, "mistuned_sol400.dat", debug=True)

# Read reference results
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

print("\nNASTRAN data import completed")
    
fig=plt.figure(figsize=(15,8))
plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[0][0][:]-1, mistuned_static_deform[0][2][:],'r',label='Gravity load' )
plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[1][0][:]-1, mistuned_static_deform[1][2][:],'b',label='10 N distributed load')
plt.plot(mistuned_grid_coords[0,:] + mistuned_static_deform[2][0][:]-1, mistuned_static_deform[2][2][:],'g',label='25 N distributed load')
plt.plot(mistuned_grid_coords[0,:]-1, mistuned_grid_coords[2,:]-1, 'k',label='No load')
plt.rcParams["lines.linewidth"] = 3
plt.title('Mistuned Static deformation with loading',fontsize=32 )
plt.xlabel("Beam length [m]",fontsize=26)
plt.ylabel("Vertical displacement [m]",fontsize=26)
plt.ax = plt.gca()
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=22)
# plt.ylim(-1,1.5)

fig.savefig("out/mistuned_beam_static_1_5m_z.svg",bbox_inches='tight')


# Load the model bdf using pynastran
# open the original .bdf file
model = BDF()
model.read_bdf("inp/refBeam.bdf", punch=True)
print(model.get_bdf_stats())
   
# get the number of properties of each type - separate mass and stiffness
elemPropKeys = list(model.properties.keys())
matStiffPropKeys = list(model.materials.keys())
matMassPropKeys = list(model.materials.keys())
conMassKeys = list(model.masses.keys())
   
numElemProps = len(elemPropKeys)
numMatStiffProps = len(matStiffPropKeys)
numMatMassProps = len(matMassPropKeys)
numConMasses = len(conMassKeys)
