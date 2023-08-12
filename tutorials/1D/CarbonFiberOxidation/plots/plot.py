from patoPlot import *

# List of data files 
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/profile/massLoss")
PATOFileList_.append("../output/"+region+"/profile/rho_s_0")
PATOFileList_.append("../output/"+region+"/profile/rho_s_1200")
PATOFileList_.append("../output/"+region+"/profile/rho_s_2400")
PATOFileList_.append("../output/"+region+"/profile/rho_s_3600")
PATOFileList_.append("../output/"+region+"/profile/Y[O2]_3600")
PATOFileList_.append("../output/"+region+"/profile/Y[N2]_3600")
PATOFileList_.append("../output/"+region+"/profile/Y[CO2]_3600")

# Read data files 
patoPlot = patoPlot(PATOFileList_)

# Set up plot parameters 
params = {
    'axes.labelsize': 14,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'text.usetex': False,
    'figure.figsize': [10, 5]
}
rcParams.update(params)

# initial thickness
initialThickness = 2.54e-2
# recession after 3600s
recession_3600s = 0.00239667

### Mass loss ###
# create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('m/m0 (-)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# set limits
xlim(0,3600)
# plot data
data=patoPlot.dataList[0]
plot(data[:,0],data[:,3],ls='solid', color = "black",linewidth=2, alpha = 0.9)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "massLoss.pdf"
savefig(outputFile_)

### Solid density ### 
# create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels
ax.set_xlabel('Depth (m)', fontsize=16)
ax.set_ylabel('Solid density (kg/m3)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# set limits
xlim(0, 2.54e-2)
# plot data
data=patoPlot.dataList[1]
plot(initialThickness-data[:,0],data[:,1],ls='solid', color = "black",linewidth=2, label = "t = 0 s", alpha = 0.9)
data=patoPlot.dataList[2]
plot(initialThickness-data[:,0],data[:,1],ls='dashed', color = "black",linewidth=2, label = "t = 1200 s", alpha = 0.9)
data=patoPlot.dataList[3]
plot(initialThickness-data[:,0],data[:,1],ls='dotted', color = "black",linewidth=2, label = "t = 2400 s", alpha = 0.9)
data=patoPlot.dataList[4]
plot(initialThickness-data[:,0],data[:,1],".-", color = "black",linewidth=2, label = "t = 3600 s", alpha = 0.9)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=4, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "solidDensity.pdf"
savefig(outputFile_)

### Mass fractions ###  
# create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels    
ax.set_xlabel('Depth (m)', fontsize=16)
ax.set_ylabel('Mass fractions (-)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# set limits
xlim(0,2.54e-2)
# plot data 
data=patoPlot.dataList[5]
plot(initialThickness-data[:,0],data[:,1],ls='solid', color = "black",linewidth=2, label = "O2", alpha = 0.9)
data=patoPlot.dataList[6]
plot(initialThickness-data[:,0],data[:,1],ls='dashed', color = "black",linewidth=2, label = "N2", alpha = 0.9)
data=patoPlot.dataList[7]
plot(initialThickness-data[:,0],data[:,1],ls='dotted', color = "black",linewidth=2, label = "CO2", alpha = 0.9)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=5, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "massFractions.pdf"
savefig(outputFile_)
