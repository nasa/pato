from patoPlot import *

# List of data files 
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/Ta_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../output/"+region+"/profile/Y[O2]_120")
PATOFileList_.append("../output/"+region+"/profile/Y[CO]_120")

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
# create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('Temperature (K)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x and y limits
ylim(200,1400)
xlim(0,120)
# plot data
data=patoPlot.dataList[0]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.9)
data=patoPlot.dataList[1]
plot(data[:,0],data[:,1],ls='solid', color = "purple",linewidth=2,label="TC1", alpha = 0.8)
plot(data[:,0],data[:,2],ls='solid', color = "red",linewidth=2,label="TC2", alpha = 0.8)
plot(data[:,0],data[:,3],ls='solid', color = "blue",linewidth=2,label="TC3", alpha = 0.8)
plot(data[:,0],data[:,4],ls='solid', color = "green",linewidth=2,label="TC4", alpha = 0.8)
plot(data[:,0],data[:,5],ls='solid', color = "grey",linewidth=2,label="TC5", alpha = 0.8)
plot(data[:,0],data[:,6],ls='solid', color = "orange",linewidth=2,label="TC6", alpha = 0.8)
plot(data[:,0],data[:,7],ls='solid', color = "yellow",linewidth=2,label="TC7", alpha = 0.8)
plt.tight_layout()
plt.grid()
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "Temperature.pdf"
savefig(outputFile_)

### Profile of mass fractions ###
# initial thickness 
initThick = 0.05
# create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels
ax.set_xlabel('Depth (m)', fontsize=16)
ax.set_ylabel('Mass fractions (-)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x and y limits
ylim(0,0.25)
xlim(0,initThick)
# plot data
data=patoPlot.dataList[2]
plot(initThick - data[:,0],data[:,1],ls='solid', color = "grey",linewidth=2,label="O2", alpha = 0.9)
data=patoPlot.dataList[3]
plot(initThick - data[:,0],data[:,1],ls='solid', color = "blue",linewidth=2,label="CO", alpha = 0.9)
plt.tight_layout()
plt.grid()
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "massFractions.pdf"
savefig(outputFile_)
