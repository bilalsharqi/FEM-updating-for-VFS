# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 20:18:09 2021

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

# Bounds function

def bounds_func(lb,ub):

    # Load the reference model bdf using pynastran
    # open the .bdf file
    ref_model = BDF()
    ref_model.read_bdf("inp/refBeam.bdf", punch=True)
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
    
    # read reference stiffness values
    ref_stiffness = np.zeros(numMatStiffProps)
    ref_mass = np.zeros(numMatMassProps)
    ref_nu = np.zeros(numMatStiffProps)
    
    for stiffness in range(len(matStiffPropKeys)):
        ref_stiffness[stiffness] = ref_model.materials[stiffness+1].e
        
    # read reference mass values
    for mass in range(len(matMassPropKeys)):
        ref_mass[mass] = ref_model.materials[mass+1].rho
        
    # read reference nu values
    for nu in range(len(matMassPropKeys)):
        ref_nu[nu] = ref_model.materials[nu+1].nu
        
    x_ref = np.array([np.transpose(ref_mass)])
    x_ref = np.append(x_ref,[np.transpose(ref_stiffness)])
    # x = np.append(x,[np.transpose(ref_nu)])
    x_lb = lb*np.float64(x_ref)
    x_ub = ub*np.float64(x_ref)
    bnds = tuple(((x_lb[i],x_ub[i]) for i in range(len(x_ref))))
    
    return bnds

# test_bnds = bounds_func(0.9, 1.1)