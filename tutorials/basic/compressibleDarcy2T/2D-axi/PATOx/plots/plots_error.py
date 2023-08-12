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
#                                          Solid Temperature                                         #
######################################################################################################

# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Set labels
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Solid Temperature error', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')

# x limit
DataTBas=data.dataList[0]
DataT1=data.dataList[1]
DataT2=data.dataList[2]
DataT3=data.dataList[3]
DataT4=data.dataList[4]
DataT5=data.dataList[5]

# Be careful, we totally do not check that the data are not 0.
DataT8= abs((DataTBas-DataT3)/DataTBas)
DataT9= abs((DataT1-DataT4)/DataT1)
DataT10= abs((DataT2-DataT5)/DataT2)

xlim = DataTBas[:,1]

# Plot
# Error at 0s is for checking that initial conditions are identical.
plot(DataTBas[:,0],DataT8[:,3],ls='solid', color = "blue",linewidth=3, alpha = 1, label='5s')
plot(DataTBas[:,0],DataT9[:,3],ls='solid', color = "green",linewidth=3, alpha = 1, label='10s')
plot(DataTBas[:,0],DataT10[:,3],ls='solid', color = "cyan",linewidth=3, alpha = 1, label='15s')

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
outputFile_ = "./plots/Ta_error.pdf"
savefig(outputFile_)


######################################################################################################
#                                          Gas Temperature                                           #
######################################################################################################

# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Set labels
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Gas Temperature error', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')

# Plot
# Error at 0s is for checking that initial conditions are identical.
plot(DataTBas[:,0],DataT8[:,2],ls='solid', color = "blue",linewidth=3, alpha = 1, label='5s')
plot(DataTBas[:,0],DataT9[:,2],ls='solid', color = "green",linewidth=3, alpha = 1, label='10s')
plot(DataTBas[:,0],DataT10[:,2],ls='solid', color = "cyan",linewidth=3, alpha = 1, label='15s')

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
outputFile_ = "./plots/Tg_error.pdf"
savefig(outputFile_)


######################################################################################################
#                                          Gas velocity                                              #
######################################################################################################

# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Set labels
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Gas velocity error', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')

# x limit
DataTBas=np.delete(data.dataList[6],[2,3],1)
DataT1=np.delete(data.dataList[7],[2,3],1)
DataT2=np.delete(data.dataList[8],[2,3],1)

DataT4=np.delete(data.dataList[9],[2,3],1)
DataT5=np.delete(data.dataList[10],[2,3],1)
DataT6=np.delete(data.dataList[11],[2,3],1)

DataT8= abs((DataTBas-DataT4)/DataTBas)
DataT9= abs((DataT1-DataT5)/DataT1)
DataT10= abs((DataT2-DataT6)/DataT2)

xlim = DataTBas[:,1]

# Plot
plot(DataTBas[:,0],DataT8[:,1],ls='solid', color = "blue",linewidth=3, alpha = 1, label='PATOx 5s')
plot(DataTBas[:,0],DataT9[:,1],ls='solid', color = "green",linewidth=3, alpha = 1, label='PATOx 10s')
plot(DataTBas[:,0],DataT10[:,1],ls='solid', color = "cyan",linewidth=3, alpha = 1, label='PATOx 15s')

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
outputFile_ = "./plots/vG_error.pdf"
savefig(outputFile_)
