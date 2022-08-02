#!/usr/bin/env python3

import patoPlot
from pylab import *

######################################################################################################
#                                            Data reading                                            #
######################################################################################################

# List of data files
PATOFileList_=[]

# Pressure/Temperatures PATOx
PATOFileList_.append("./postProcessing/singleGraph/porousMat/5/line_p_Tg_Ta.xy")
PATOFileList_.append("./postProcessing/singleGraph/porousMat/10/line_p_Tg_Ta.xy")
PATOFileList_.append("./postProcessing/singleGraph/porousMat/15/line_p_Tg_Ta.xy")

# Pressure/Temperatures compressibleDarcy2T
PATOFileList_.append("../postProcessing/singleGraph/5/line_p_T_g_T_s.xy")
PATOFileList_.append("../postProcessing/singleGraph/10/line_p_T_g_T_s.xy")
PATOFileList_.append("../postProcessing/singleGraph/15/line_p_T_g_T_s.xy")

# Gas velocity PATOx
PATOFileList_.append("./postProcessing/singleGraph/porousMat/5/line_vG.xy")
PATOFileList_.append("./postProcessing/singleGraph/porousMat/10/line_vG.xy")
PATOFileList_.append("./postProcessing/singleGraph/porousMat/15/line_vG.xy")

# Gas velocity compressibleDarcy2T
PATOFileList_.append("../postProcessing/singleGraph/5/line_u_g.xy")
PATOFileList_.append("../postProcessing/singleGraph/10/line_u_g.xy")
PATOFileList_.append("../postProcessing/singleGraph/15/line_u_g.xy")

# Read data files
data = patoPlot.patoPlot(PATOFileList_)


######################################################################################################
#                                           Set up plot parameters                                   #
######################################################################################################

params = {
    'axes.labelsize': 14,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'text.usetex': False,
    'figure.figsize': [10, 5]
}
rcParams.update(params)


######################################################################################################
#                                          Pressure                                                  #
######################################################################################################

# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Set labels
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Pressure (Pa)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')

# x limit
DataTBas=data.dataList[0]
DataT1=data.dataList[1]
DataT2=data.dataList[2]
DataT3=data.dataList[3]
DataT4=data.dataList[4]
DataT5=data.dataList[5]

xlim = DataTBas[:,1]

# Plot
plot(DataTBas[:,0],DataTBas[:,1],ls='solid', color = "blue",linewidth=3, alpha = 1, label='PATOx 5s')
plot(DataTBas[:,0],DataT1[:,1],ls='solid', color = "green",linewidth=3, alpha = 1, label='PATOx 10s')
plot(DataTBas[:,0],DataT2[:,1],ls='solid', color = "cyan",linewidth=3, alpha = 1, label='PATOx 15s')
plot(DataTBas[::3,0],DataT3[::3,1],'ko',color = "blue", markersize=4, marker="v", mfc='none',label='5s')
plot(DataTBas[::3,0],DataT4[::3,1],'ko', color = "green", markersize=4, marker="^", mfc='none',label='10s')
plot(DataTBas[::3,0],DataT5[::3,1],'ko', color = "cyan", markersize=4, marker="s", mfc='none',label='15s')

# Tight layout
plt.tight_layout()

# Add grid
plt.grid()

# Put legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')

# Write figure
outputFile_ = "./plots/p.pdf"
savefig(outputFile_)


######################################################################################################
#                                          Solid Temperature                                         #
######################################################################################################

# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Set labels
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Solid Temperature (K)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')

# Plot
plot(DataTBas[:,0],DataTBas[:,3],ls='solid', color = "blue",linewidth=3, alpha = 1, label='PATOx 5s')
plot(DataTBas[:,0],DataT1[:,3],ls='solid', color = "green",linewidth=3, alpha = 1, label='PATOx 10s')
plot(DataTBas[:,0],DataT2[:,3],ls='solid', color = "cyan",linewidth=3, alpha = 1, label='PATOx 15s')
plot(DataTBas[::3,0],DataT3[::3,3],'ko',color = "blue", markersize=4, marker="v", mfc='none',label='5s')
plot(DataTBas[::3,0],DataT4[::3,3],'ko', color = "green", markersize=4, marker="^", mfc='none',label='10s')
plot(DataTBas[::3,0],DataT5[::3,3],'ko', color = "cyan", markersize=4, marker="s", mfc='none',label='15s')

# Tight layout
plt.tight_layout()

# Add grid
plt.grid()

# Put legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')

# Write figure
outputFile_ = "./plots/Ta.pdf"
savefig(outputFile_)


######################################################################################################
#                                          Gas Temperature                                           #
######################################################################################################

# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Set labels
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Gas Temperature (K)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')

# Plot
plot(DataTBas[:,0],DataTBas[:,2],ls='solid', color = "blue",linewidth=3, alpha = 1, label='PATOx 5s')
plot(DataTBas[:,0],DataT1[:,2],ls='solid', color = "green",linewidth=3, alpha = 1, label='PATOx 10s')
plot(DataTBas[:,0],DataT2[:,2],ls='solid', color = "cyan",linewidth=3, alpha = 1, label='PATOx 15s')
plot(DataTBas[::3,0],DataT3[::3,2],'ko',color = "blue", markersize=4, marker="v", mfc='none',label='5s')
plot(DataTBas[::3,0],DataT4[::3,2],'ko', color = "green", markersize=4, marker="^", mfc='none',label='10s')
plot(DataTBas[::3,0],DataT5[::3,2],'ko', color = "cyan", markersize=4, marker="s", mfc='none',label='15s')

# Tight layout
plt.tight_layout()

# Add grid
plt.grid()

# Put legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')

# Write figure
outputFile_ = "./plots/Tg.pdf"
savefig(outputFile_)


######################################################################################################
#                                          Gas velocity                                              #
######################################################################################################

# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Set labels
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Gas velocity (m/s)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')

# x limit
DataTBas=data.dataList[6]
DataT1=data.dataList[7]
DataT2=data.dataList[8]
DataT3=data.dataList[9]
DataT4=data.dataList[10]
DataT5=data.dataList[11]

xlim = DataTBas[:,1]

# Plot
plot(DataTBas[:,0],DataTBas[:,1],ls='solid', color = "blue",linewidth=3, alpha = 1, label='PATOx 5s')
plot(DataTBas[:,0],DataT1[:,1],ls='solid', color = "green",linewidth=3, alpha = 1, label='PATOx 10s')
plot(DataTBas[:,0],DataT2[:,1],ls='solid', color = "cyan",linewidth=3, alpha = 1, label='PATOx 15s')
plot(DataTBas[::3,0],DataT3[::3,1],'ko',color = "blue", markersize=4, marker="v", mfc='none',label='5s')
plot(DataTBas[::3,0],DataT4[::3,1],'ko', color = "green", markersize=4, marker="^", mfc='none',label='10s')
plot(DataTBas[::3,0],DataT5[::3,1],'ko', color = "cyan", markersize=4, marker="s", mfc='none',label='15s')

# Tight layout
plt.tight_layout()

# Add grid
plt.grid()

# Put legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')

# Write figure
outputFile_ = "./plots/vG.pdf"
savefig(outputFile_)
