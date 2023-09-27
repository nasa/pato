from patoPlot import *

# List of data files
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/Ta_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../data/ref/FIAT/T")
PATOFileList_.append("../output/"+region+"/mass")
PATOFileList_.append("../data/ref/FIAT/pyrolysisFront")
PATOFileList_.append("../output/"+region+"/massLoss")

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
data=patoPlot.dataList[2]
plot(data[:,0],data[:,1],ls='dashed', color = "black",linewidth=2, alpha = 0.8)
plot(data[:,0],data[:,2],ls='dashed', color = "purple",linewidth=2, alpha = 0.8)
plot(data[:,0],data[:,3],ls='dashed', color = "red",linewidth=2, alpha = 0.8)
plot(data[:,0],data[:,4],ls='dashed', color = "blue",linewidth=2, alpha = 0.8)
plot(data[:,0],data[:,5],ls='dashed', color = "green",linewidth=2, alpha = 0.8)
plot(data[:,0],data[:,6],ls='dashed', color = "orange",linewidth=2, alpha = 0.8)
# plot the text
plot([81, 86], [1730, 1730], lw=2, ls="solid", color="black")
plot([81, 86], [1645, 1645], lw=2, ls="dashed", color="black")
text(88, 1705, 'PATO', fontsize=12)
text(88, 1620, 'FIAT', fontsize=12)
plt.tight_layout()
plt.grid()
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "Temperature.pdf"
savefig(outputFile_)

### Pyrolysis mass flux ### 
# Create figure 
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels 
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('$\dot{m}_g$ (kg/m2/s)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,60)
ylim(0,0.06)
# plot the data 
data=patoPlot.dataList[3]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="", alpha = 0.8)
data=patoPlot.dataList[4]
plot(data[:,0],data[:,2],ls='dashed', color = "black",linewidth=2,label="", alpha = 0.8)
# plot the text           
plot([31, 33], [0.056, 0.056], lw=2, ls="solid", color="black")
plot([31, 33], [0.052, 0.052], lw=2, ls="dashed", color="black")
text(35, 0.055, 'PATO', fontsize=12)
text(35, 0.051, 'FIAT', fontsize=12)
plt.tight_layout()
plt.grid()
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "PyrolysisMassFlux.pdf"
savefig(outputFile_)

### 98% virgin / 2% char ### 
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels    
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('Depth (m)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit 
xlim(0,60)
ylim(0,0.025)
# plot the data 
data=patoPlot.dataList[3]
plot(data[:,0],data[:,3],ls='solid', color = "black",linewidth=2,label="98% virgin", alpha = 0.8)
plot(data[:,0],data[:,4],ls='solid', color = "grey",linewidth=2,label="2% char", alpha = 0.8)
data=patoPlot.dataList[4]
plot(data[:,0],data[:,8],ls='dashed', color = "black",linewidth=2,label="", alpha = 0.8)
plot(data[:,0],data[:,7],ls='dashed', color = "grey",linewidth=2,label="", alpha = 0.8)
# plot the text  
plot([31, 33], [0.0235, 0.0235], lw=2, ls="solid", color="black")
plot([31, 33], [0.0215, 0.0215], lw=2, ls="dashed", color="black")
text(35, 0.023, 'PATO', fontsize=12)
text(35, 0.021, 'FIAT', fontsize=12)
plt.tight_layout()
plt.grid()
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "virgin_char_fronts.pdf"
savefig(outputFile_)
