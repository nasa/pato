from patoPlot import *

# List of data files
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/Ta_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../output/porousMat/profile/eps_s_phase1_120")
PATOFileList_.append("../output/porousMat/profile/eps_s_phase2_120")
PATOFileList_.append("../output/porousMat/profile/eps_s_phase3_120")

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

### Temperature ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('Temperature (K)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
# plot the data
data=patoPlot.dataList[0]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
data=patoPlot.dataList[1]
plot(data[:,0],data[:,1],ls='solid', color = "purple",linewidth=2,label="TC1", alpha = 0.8)
plot(data[:,0],data[:,2],ls='solid', color = "red",linewidth=2,label="TC2", alpha = 0.8)
plot(data[:,0],data[:,3],ls='solid', color = "blue",linewidth=2,label="TC3", alpha = 0.8)
plot(data[:,0],data[:,4],ls='solid', color = "green",linewidth=2,label="TC4", alpha = 0.8)
plot(data[:,0],data[:,5],ls='solid', color = "grey",linewidth=2,label="TC5", alpha = 0.8)
plot(data[:,0],data[:,6],ls='solid', color = "orange",linewidth=2,label="TC6", alpha = 0.8)
plot(data[:,0],data[:,7],ls='solid', color = "yellow",linewidth=2,label="TC7", alpha = 0.8)
# tight layout
plt.tight_layout()
# add grid
plt.grid()
# put legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "Temperature.pdf"
savefig(outputFile_)

### Solid volume fraction ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Depth (cm)', fontsize=16)
ax.set_ylabel('Solid volume fraction (-)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x and y limits 
ylim(0,0.11)
xlim(0,2)
# plot the data
data=patoPlot.dataList[2]
plot(100*(0.05-data[:,0]),data[:,1],ls='dashed', label="phase 1 (fiber) and phase 2 (matrix)", color = "black",linewidth=2, alpha = 0.9)
data=patoPlot.dataList[4]
plot(100*(0.05-data[:,0]),data[:,1],ls='solid',label="phase 3 (grading)", color = "black",linewidth=2, alpha = 0.9)
# tight layout
plt.tight_layout()
# add grid
plt.grid()
# put legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=5, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "SolidVolumeFraction.pdf"
savefig(outputFile_)

