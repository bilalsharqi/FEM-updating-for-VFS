# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 16:37:44 2023

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
# from pyMSCNastranUtilities import *
plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

# frequency and error comparisons

Errors = np.array([[11.37,0.47,4.52,0.12,0.16],
                    [11.70,0.04,0.66,0.25,0.03],
                    [10.22,0.09,1.75,0.88,0.27],
                    [12.10,1.91,3.66,0.89,1.31],
                    [7.12,1.04,5.99,1.04,1.34]])


# bar plot for errors
labels = ['1', '2', '3', '4', '5']
lab_gvt_ini = Errors[0:,0]
lab_gvt_fin = Errors[0:,1]
wt_gvt_ini = Errors[0:,2]
wt_gvt_fin = Errors[0:,3]
lab_gvt_fin_split = Errors[0:,4]

# figure 1 Lab GVT
x = np.arange(len(labels))+1  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(figsize=(16,9))

rects1 = ax.bar(x, lab_gvt_ini, width, label='Shaker table GVT initial',alpha=0.4)
rects2 = ax.bar(x, lab_gvt_fin, width, label='Shaker table GVT final')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Error [%]',fontsize=26)
ax.set_xlabel('Mode number',fontsize=26)
# ax.set_title('Errors')
ax.set_xticks(x)
ax.legend(loc="upper left")

plt.ax = plt.gca()
plt.xticks(fontsize=25)

plt.yticks(fontsize=25)
plt.legend(fontsize=22, frameon=False)
plt.ax.spines['top'].set_visible(False)
plt.ax.spines['right'].set_visible(False)
fig.savefig('mode_shapes/Errors_ini_fin_lab.svg',bbox_inches='tight')

fig.tight_layout()

plt.show()

# figure 2 - WT GVT
x = np.arange(len(labels))+1  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(figsize=(16,9))

rects1 = ax.bar(x, wt_gvt_ini, width, label='Wind tunnel GVT initial',alpha=0.4)
rects2 = ax.bar(x, wt_gvt_fin, width, label='Wind tunnel GVT final')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Error [%]',fontsize=26)
ax.set_xlabel('Mode number',fontsize=26)
# ax.set_title('Errors')
ax.set_xticks(x)
ax.legend(loc="upper left")

plt.ax = plt.gca()
plt.xticks(fontsize=25)

plt.yticks(fontsize=25)
plt.legend(fontsize=22, frameon=False)
plt.ax.spines['top'].set_visible(False)
plt.ax.spines['right'].set_visible(False)
fig.savefig('mode_shapes/Errors_ini_fin_wtt.svg',bbox_inches='tight')

fig.tight_layout()

plt.show()

# figure 3 - Lab GVT Split Segments
x = np.arange(len(labels))+1  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(figsize=(16,9))

rects1 = ax.bar(x, lab_gvt_ini, width, label='Initial FEM',alpha=0.4)
rects2 = ax.bar(x, lab_gvt_fin, width, label='Updated FEM - constant properties')
rects3 = ax.bar(x, lab_gvt_fin_split, width, label='Updated FEM - variable properties')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Error [%]',fontsize=26)
ax.set_xlabel('Mode number',fontsize=26)
# ax.set_title('Errors')
ax.set_xticks(x)

# ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

plt.ax = plt.gca()
plt.xticks(fontsize=25)
plt.ax.legend(loc="best")
plt.yticks(fontsize=25)
plt.legend(fontsize=22, frameon=False)
plt.ax.spines['top'].set_visible(False)
plt.ax.spines['right'].set_visible(False)
fig.savefig('mode_shapes/Errors_ini_fin_lab_split.svg',bbox_inches='tight')

fig.tight_layout()

plt.show()

