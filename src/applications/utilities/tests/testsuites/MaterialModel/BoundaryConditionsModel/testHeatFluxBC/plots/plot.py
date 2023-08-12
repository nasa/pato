from patoPlot import *
matplotlib.use('agg')

# List of data files
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/mDotGw_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../output/"+region+"/scalar/Ta_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/h_r_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/rhoeUeCH_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/eps_g_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/mu_g_surfacePatch")
PATOFileList_.append("../output/"+region+"/tensor/Kyy_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/mDotCw_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/p_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/p_plot")
PATOFileList_.append("../output/"+region+"/scalar/rho_g_surfacePatch")

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

### mDotGw ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('mDotGw', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[0]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "mDotGw.pdf"
savefig(outputFile_)

### Temperature ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('T', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[2]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
data=patoPlot.dataList[1]
plot(data[:,0],data[:,1],ls='solid', color = "purple",linewidth=2,label="TC1", alpha = 0.8)
plot(data[:,0],data[:,2],ls='solid', color = "red",linewidth=2,label="TC2", alpha = 0.8)
plot(data[:,0],data[:,3],ls='solid', color = "blue",linewidth=2,label="TC3", alpha = 0.8)
plot(data[:,0],data[:,4],ls='solid', color = "green",linewidth=2,label="TC4", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "Temperature.pdf"
savefig(outputFile_)

### h_r ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('h_r', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[3]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "h_r.pdf"
savefig(outputFile_)

### Ch ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('Ch', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[4]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "Ch.pdf"
savefig(outputFile_)

### eps_g ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('eps_g', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[5]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "eps_g.pdf"
savefig(outputFile_)

### mu_g ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('mu_g', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[6]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "mu_g.pdf"
savefig(outputFile_)

### K ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('K', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[7]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "K.pdf"
savefig(outputFile_)

### mDotCw ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('mDotCw', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[8]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "mDotCw.pdf"
savefig(outputFile_)

### p ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('p', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[9]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
data=patoPlot.dataList[10]
plot(data[:,0],data[:,1],ls='solid', color = "purple",linewidth=2,label="First Cell", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "p.pdf"
savefig(outputFile_)

### dp ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('dp', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[9]
for i in range(len(data[:,0])):
    data[i,1]=data[i,1]-patoPlot.dataList[10][i,1]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "dp.pdf"
savefig(outputFile_)


### rho_g ###
# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Time (s)', fontsize=16)
ax.set_ylabel('rho_g', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0,120)
#ylim(200,1700)
# plot the data
data=patoPlot.dataList[11]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="Wall", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.tight_layout()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "rho_g.pdf"
savefig(outputFile_)
