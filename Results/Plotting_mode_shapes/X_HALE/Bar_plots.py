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

Errors = np.array([[4.35,	4.17,	1.65,	1.65],
                    [3.46,	0.06,	3.96,	0.11],
                    [0.57,	1.95,	2.93,	0.25],
                    [1.33,	0.84,	0.39,	0.15],
                    [1.49,	1.54,	1.83,	1.40],
                    [1.24,	0.65,	1.84,	2.01],
                    [8.21,	0.22,	5.97,	1.82],
                    [7.81,	1.98,	4.24,	1.51]])


# bar plot for errors
labels = ['1', '2', '3', '4', '5','6', '7', '8']
inboard_ini = Errors[0:,0]
inboard_fin = Errors[0:,1]
outboard_ini = Errors[0:,2]
outboard_fin = Errors[0:,3]

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(figsize=(16,9))

rects1 = ax.bar(x - width/2, np.abs(inboard_fin - inboard_ini), width, label='Inboard')
rects2 = ax.bar(x + width/2, np.abs(outboard_fin - outboard_ini), width, label='Outboard')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Reduction in percent error',fontsize=26)
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
fig.savefig('mode_shapes/Reduction_Errors_ini_fin_XHALE.svg',bbox_inches='tight')

fig.tight_layout()

plt.show()

# absolute errors
# bar plot for errors
labels = [ '1','2', '3', '4', '5','6', '7', '8']
inboard_ini = Errors[0:,0]
inboard_fin = Errors[0:,1]
outboard_ini = Errors[0:,2]
outboard_fin = Errors[0:,3]

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(figsize=(16,9))

rects1 = ax.bar(x - width/2, inboard_ini, width, label='Inboard initial', alpha=0.4)
rects2 = ax.bar(x + width/2, outboard_ini, width, label='Outboard initial', alpha=0.4)

rects3 = ax.bar(x - width/2, inboard_fin, width, label='Inboard final')
rects4 = ax.bar(x + width/2, outboard_fin, width, label='Outboard final')

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
fig.savefig('mode_shapes/Errors_ini_fin_XHALE.svg',bbox_inches='tight')

fig.tight_layout()

plt.show()

# # models calibrated using complex load linear vs nonlinear optimization
# lin_v_nl = np.array([[1,	0.68,	0.68,	0.83,	0.68,	0.23,],
#         [2,	2.87,	2.87,	0.03,	2.88,	0.10,],
#         [3,	4.25,	4.27,	0.42,	4.24,	0.21,],
#         [4,	11.71,	11.73,	0.24,	11.69,	0.16,],
#         [5,	17.20,	17.20,	0.04,	17.20,	0.04,],
#         [6,	18.04,	18.03,	0.05,	18.05,	0.02,],
#         [7,	22.18,	22.22,	0.17,	22.16,	0.11,],
#         [8,	34.81,	34.91,	0.26,	34.80,	0.04,],
#         [9,	47.90,	47.87,	0.05,	47.93,	0.07,],
#         [10,	48.23,	48.58,	0.72,	48.25,	0.05,],
#         [11,	52.23,	52.24,	0.01,	52.18,	0.09,],
#         [12,	54.89,	58.01,	5.69,	55.07,	0.33,],
#         [13,	60.57,	62.54,	3.26,	60.64,	0.12,],
#         [14,	63.92,	64.00,	0.12,	63.88,	0.06,],
#         [15,	67.16,	69.97,	4.19,	67.25,	0.14]])

# # bar plot for errors
# labels = ['1', '2', '3', '4', '5','6', '7', '8', '9', '10','11', '12', '13', '14', '15']
# lin_err = lin_v_nl[0:,3]
# nl_err = lin_v_nl[0:,5]


# x = np.arange(len(labels))+1  # the label locations
# width = 0.35  # the width of the bars

# fig, ax = plt.subplots(figsize=(16,9))
# # fig=plt.figure(figsize=(16,9))
# # rects1 = ax.bar(x - width, no_load_ini, width, label='No load initial')
# # rects2 = ax.bar(x, g_ini, width, label='- 1 g load initial')
# # rects3 = ax.bar(x + width, g_2_5_ini, width, label='2.5 g load initial')

# rects1 = ax.bar(x - width/2, lin_err, width, label='Conventional optimization vs. reference FEM')
# rects2 = ax.bar(x + width/2, nl_err, width, label='Modified optimzation vs. reference FEM')

# # rects4 = ax.bar(x + 2*width, no_load_fin, width, label='No load final')
# # rects5 = ax.bar(x + 3*width, g_fin, width, label='- 1 g load final')
# # rects6 = ax.bar(x + 4*width, g_2_5_fin, width, label='2.5 g load final')

# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Error [%]',fontsize=26)
# ax.set_xlabel('Mode number',fontsize=26)
# # ax.set_title('Errors')
# ax.set_xticks(x)
# ax.legend()

# # plt.xlabel('Mode number',fontsize=26)
# # plt.rcParams["lines.linewidth"] = 3
# # plt.ylabel('Error [%]',fontsize=26)
# plt.ax = plt.gca()
# plt.xticks(fontsize=25)
# # plt.ylim(0,72)
# # plt.xticks(np.arange(1, n_freq+1, step=1.0))
# # plt.xticks([r + barWidth for r in range(len(data[:,1]))],
#     # ['1', '2', '3', '4', '5','6', '7', '8', '9', '10','11', '12', '13', '14', '15'])
# plt.yticks(fontsize=25)
# plt.legend(fontsize=22, frameon=False)
# plt.ax.spines['top'].set_visible(False)
# plt.ax.spines['right'].set_visible(False)
# fig.savefig('mode_shapes/linear_v_nonlinear_modes_bar.svg',bbox_inches='tight')

# fig.tight_layout()

# plt.show()

# # plot for beam errors
# beam_errors = np.array([[5.16,	3.18,	0.85,	5.67,	0.14,	0.15,	0.16,	0.14],
#                        [6.23,	14.9,	8.88,	6.14,	0.03,	0.07,	0.04,	0.04],
#                        [12.8,	6.44,	5.81,	4.30,	0.05,	0.00,	0.06,	0.01],
#                        [5.58,	6.26,	6.28,	5.36,	0.07,	0.06,	0.02,	0.07],
#                        [5.55,	5.96,	6.08,	5.41,	0.05,	0.02,	0.03,	0.07],
#                        [5.76,	10.8,	10.3,	3.54,	0.01,	0.10,	0.15,	0.03],
#                        [4.58,	6.03,	6.09,	5.66,	0.01,	0.02,	0.04,	0.01],
#                        [5.50,	5.69,	5.74,	5.43,	0.05,	0.05,	0.05,	0.05],
#                        [3.14,	5.67,	5.74,	4.74,	0.14,	0.06,	0.07,	0.01],
#                        [5.51,	0.60,	3.93,	5.46,	0.06,	0.17,	0.09,	0.06],
#                        [5.64,	5.76,	5.80,	5.60,	0.03,	0.03,	0.03,	0.03],
#                        [4.39,	5.57,	5.60,	3.15,	0.05,	0.05,	0.04,	0.20],
#                        [1.09,	4.73,	5.61,	5.45,	0.04,	0.36,	0.04,	0.05],
#                        [5.50,	5.58,	3.80,	5.47,	0.05,	0.05,	0.26,	0.05],
#                        [4.08,	1.97,	5.67,	3.07,	0.17,	0.03,	0.02,	0.24]])

# # bar plot for errors
# labels = ['1', '2', '3', '4', '5','6', '7', '8', '9', '10','11', '12', '13', '14', '15']
# grav_ini = beam_errors[0:,0]
# grav_fin = beam_errors[0:,4]
# ten_ini = beam_errors[0:,1]
# ten_fin = beam_errors[0:,5]
# twenty_five_ini = beam_errors[0:,2]
# twenty_five_fin = beam_errors[0:,6]
# no_load_ini = beam_errors[0:,3]
# no_load_fin = beam_errors[0:,7]

# x = np.arange(len(labels))+1  # the label locations
# width = 0.15  # the width of the bars

# fig, ax = plt.subplots(figsize=(16,9))

# rects4 = ax.bar(x + 3*width/2, np.abs(no_load_fin - no_load_ini), width, label='No load')
# rects1 = ax.bar(x - 3*width/2, np.abs(grav_fin - grav_ini), width, label='Self-weight')
# rects2 = ax.bar(x - width/2, np.abs(ten_fin - ten_ini), width, label='10-N distributed load')
# rects3 = ax.bar(x + width/2, np.abs(twenty_five_fin - twenty_five_ini), width, label='25-N distributed load')


# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Reduction in error [%]',fontsize=26)
# ax.set_xlabel('Mode number',fontsize=26)
# # ax.set_title('Errors')
# ax.set_xticks(x)
# ax.legend()

# # plt.xlabel('Mode number',fontsize=26)
# # plt.rcParams["lines.linewidth"] = 3
# # plt.ylabel('Error [%]',fontsize=26)
# plt.ax = plt.gca()
# plt.xticks(fontsize=25)
# # plt.ylim(0,72)
# # plt.xticks(np.arange(1, n_freq+1, step=1.0))
# # plt.xticks([r + barWidth for r in range(len(data[:,1]))],
#     # ['1', '2', '3', '4', '5','6', '7', '8', '9', '10','11', '12', '13', '14', '15'])
# plt.yticks(fontsize=25)
# plt.legend(fontsize=22, frameon=False)
# plt.ax.spines['top'].set_visible(False)
# plt.ax.spines['right'].set_visible(False)
# fig.savefig('mode_shapes/Error_cases_all_beam.svg',bbox_inches='tight')

# fig.tight_layout()

# plt.show()

# # plot for beam errors with static deflection in objective function
# # plot for beam errors
# beam_errors = np.array([[5.16,	3.18,	0.85,	5.67,	0.04,	0.10,	0.14,	0.06],
#                         [6.23,	14.9,	8.88,	6.14,	0.31,	3.01,	2.94,	0.29],
#                         [12.8,	6.44,	5.81,	4.30,	0.34,	0.34,	0.40,	4.98],
#                         [5.58,	6.26,	6.28,	5.36,	0.34,	0.35,	0.35,	5.18],
#                         [5.55,	5.96,	6.08,	5.41,	0.33,	0.44,	0.58,	0.30],
#                         [5.76,	10.8,	10.3,	3.54,	0.12,	0.41,	1.40,	0.23],
#                         [4.58,	6.03,	6.09,	5.66,	4.63,	0.18,	0.23,	0.10],
#                         [5.50,	5.69,	5.74,	5.43,	0.29,	0.31,	0.31,	0.29],
#                         [3.14,	5.67,	5.74,	4.74,	2.10,	0.27,	0.28,	0.36],
#                         [5.51,	0.60,	3.93,	5.46,	0.26,	0.84,	0.68,	0.26],
#                         [5.64,	5.76,	5.80,	5.60,	0.15,	0.15,	0.15,	0.15],
#                         [4.39,	5.57,	5.60,	3.15,	0.27,	0.28,	0.28,	2.10],
#                         [1.09,	4.73,	5.61,	5.45,	0.12,	2.56,	0.28,	0.27],
#                         [5.50,	5.58,	3.80,	5.47,	0.25,	0.27,	1.78,	0.25],
#                         [4.08,	1.97,	5.67,	3.07,	3.56,	0.19,	0.19,	1.78]])

# # bar plot for errors
# labels = ['1', '2', '3', '4', '5','6', '7', '8', '9', '10','11', '12', '13', '14', '15']
# grav_ini = beam_errors[0:,0]
# grav_fin = beam_errors[0:,4]
# ten_ini = beam_errors[0:,1]
# ten_fin = beam_errors[0:,5]
# twenty_five_ini = beam_errors[0:,2]
# twenty_five_fin = beam_errors[0:,6]
# no_load_ini = beam_errors[0:,3]
# no_load_fin = beam_errors[0:,7]

# x = np.arange(len(labels))+1  # the label locations
# width = 0.15  # the width of the bars

# fig, ax = plt.subplots(figsize=(16,9))

# rects4 = ax.bar(x + 3*width/2, np.abs(no_load_fin - no_load_ini), width, label='No load')
# rects1 = ax.bar(x - 3*width/2, np.abs(grav_fin - grav_ini), width, label='Self-weight')
# rects2 = ax.bar(x - width/2, np.abs(ten_fin - ten_ini), width, label='10-N distributed load')
# rects3 = ax.bar(x + width/2, np.abs(twenty_five_fin - twenty_five_ini), width, label='25-N distributed load')


# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Reduction in error [%]',fontsize=26)
# ax.set_xlabel('Mode number',fontsize=26)
# # ax.set_title('Errors')
# ax.set_xticks(x)
# ax.legend()

# # plt.xlabel('Mode number',fontsize=26)
# # plt.rcParams["lines.linewidth"] = 3
# # plt.ylabel('Error [%]',fontsize=26)
# plt.ax = plt.gca()
# plt.xticks(fontsize=25)
# # plt.ylim(0,72)
# # plt.xticks(np.arange(1, n_freq+1, step=1.0))
# # plt.xticks([r + barWidth for r in range(len(data[:,1]))],
#     # ['1', '2', '3', '4', '5','6', '7', '8', '9', '10','11', '12', '13', '14', '15'])
# plt.yticks(fontsize=25)
# plt.legend(fontsize=22, frameon=False)
# plt.ax.spines['top'].set_visible(False)
# plt.ax.spines['right'].set_visible(False)
# fig.savefig('mode_shapes/Error_reduction_static_def_beam.svg',bbox_inches='tight')

# fig.tight_layout()

# plt.show()

# # model calibrated using linear, unloaded and undeformed optimization
# data_103 = np.array([[0.720,	0.740,	0.720,	3.55,	0.13,],
#                      [2.900,	3.000,	2.900,	3.38,	0.13,],
#                      [4.410,	4.570,	4.410,	3.54,	0.13,],
#                      [12.02,	12.45,	12.01,	3.53,	0.12,],
#                      [17.52,	18.14,	17.52,	3.52,	0.05,],
#                      [17.78,	18.38,	17.75,	3.38,	0.13,],
#                      [22.66,	23.46,	22.64,	3.52,	0.11,],
#                      [35.65,	36.91,	35.62,	3.54,	0.08,],
#                      [48.18,	49.82,	48.12,	3.39,	0.13,],
#                      [50.13,	51.93,	50.12,	3.60,	0.02,],
#                      [52.18,	54.05,	52.18,	3.57,	0.01,],
#                      [65.20,	67.61,	65.24,	3.70,	0.06,],
#                      [66.84,	69.18,	66.76,	3.51,	0.11,],
#                      [80.01,	83.11,	80.15,	3.86,	0.17,],
#                      [85.56,	88.72,	85.61,	3.70,	0.07,]])

# # bar plot for errors
# labels = ['1', '2', '3', '4', '5','6', '7', '8', '9', '10','11', '12', '13', '14', '15']
# ini_err = data_103[0:,3]
# fin_err = data_103[0:,4]


# x = np.arange(len(labels))+1  # the label locations
# width = 0.35  # the width of the bars

# fig, ax = plt.subplots(figsize=(16,9))
# # fig=plt.figure(figsize=(16,9))
# # rects1 = ax.bar(x - width, no_load_ini, width, label='No load initial')
# # rects2 = ax.bar(x, g_ini, width, label='- 1 g load initial')
# # rects3 = ax.bar(x + width, g_2_5_ini, width, label='2.5 g load initial')

# rects1 = ax.bar(x - width/2, ini_err, width, label='Initial vs. reference FEM')
# rects2 = ax.bar(x + width/2, fin_err, width, label='Calibrated vs. reference FEM')

# # rects4 = ax.bar(x + 2*width, no_load_fin, width, label='No load final')
# # rects5 = ax.bar(x + 3*width, g_fin, width, label='- 1 g load final')
# # rects6 = ax.bar(x + 4*width, g_2_5_fin, width, label='2.5 g load final')

# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Error [%]',fontsize=26)
# ax.set_xlabel('Mode number',fontsize=26)
# # ax.set_title('Errors')
# ax.set_xticks(x)
# ax.legend()

# # plt.xlabel('Mode number',fontsize=26)
# # plt.rcParams["lines.linewidth"] = 3
# # plt.ylabel('Error [%]',fontsize=26)
# plt.ax = plt.gca()
# plt.xticks(fontsize=25)
# plt.ylim(0,5)
# # plt.xticks(np.arange(1, n_freq+1, step=1.0))
# # plt.xticks([r + barWidth for r in range(len(data[:,1]))],
#     # ['1', '2', '3', '4', '5','6', '7', '8', '9', '10','11', '12', '13', '14', '15'])
# plt.yticks(fontsize=25)
# plt.legend(fontsize=22, frameon=False)
# plt.ax.spines['top'].set_visible(False)
# plt.ax.spines['right'].set_visible(False)
# fig.savefig('mode_shapes/linear_ini_fin.svg',bbox_inches='tight')

# fig.tight_layout()

# plt.show()