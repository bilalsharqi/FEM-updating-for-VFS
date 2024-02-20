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

def bounds_func(x_ref,lb,ub):

    # apply bounds and format as tuple of tuples
    x_lb = lb*np.float64(x_ref)
    x_ub = ub*np.float64(x_ref)
    bnds = tuple(((x_lb[i],x_ub[i]) for i in range(len(x_ref))))
    
    return bnds, x_lb, x_ub

def bounds_func_single(x_ref,lb,ub):

    # apply bounds and format as tuple of tuples
    x_lb = lb*np.float64(x_ref)
    x_ub = ub*np.float64(x_ref)
    bnds = tuple(((x_lb,x_ub)))
    
    return bnds, x_lb, x_ub