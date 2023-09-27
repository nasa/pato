from patoPlot import *

# List of data files
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../plots/probe1")
PATOFileList_.append("../plots/probe2")
PATOFileList_.append("../plots/probe3")
PATOFileList_.append("../plots/surface")
PATOFileList_.append("../output/"+region+"/scalar/Ta_surfacePatch")
PATOFileList_.append("../output/"+region+"/scalar/Tg_plot")
PATOFileList_.append("../output/"+region+"/scalar/Tg_surfacePatch")

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
xlim(0,10)
ylim(273,2000)
# plot the data
# experiment
datasurface=patoPlot.dataList[4]
plot(datasurface[:,0],datasurface[:,1],'o', mfc="none", color = "black",linewidth=1,label="TCsurface_exp", alpha = 0.8)
dataprobe1=patoPlot.dataList[1]
plot(dataprobe1[:,0],dataprobe1[:,1],'o', mfc="none", color = "purple",linewidth=1,label="TC1_exp", alpha = 0.8)
dataprobe2=patoPlot.dataList[2]
plot(dataprobe2[:,0],dataprobe2[:,1],'o', mfc="none", color = "red",linewidth=1,label="TC2_exp", alpha = 0.8)
dataprobe3=patoPlot.dataList[3]
plot(dataprobe3[:,0],dataprobe3[:,1],'o', mfc="none", color = "blue",linewidth=1,label="TC3_exp", alpha = 0.8)
# simulation
dataTg=patoPlot.dataList[6]
datasurfaceg=patoPlot.dataList[7]
plot(datasurfaceg[:,0],datasurfaceg[:,1],ls='dashed', color = "black",linewidth=1,label="TCsurface_g_PATO", alpha = 0.8)
plot(dataTg[:,0],dataTg[:,1],ls='dashed', color = "purple",linewidth=1,label="TC1_g_PATO", alpha = 0.8)
plot(dataTg[:,0],dataTg[:,2],ls='dashed', color = "red",linewidth=1,label="TC2_g_PATO", alpha = 0.8)
plot(dataTg[:,0],dataTg[:,3],ls='dashed', color = "blue",linewidth=1,label="TC3_g_PATO", alpha = 0.8)
#plot(dataTg[:,0],dataTg[:,4],ls='dashed', color = "yellow",linewidth=1,label="TCbottom_g", alpha = 0.8)
dataTa=patoPlot.dataList[0]
plot(dataTa[:,0],dataTa[:,1],ls='solid', color = "purple",linewidth=1,label="TC1_s_PATO", alpha = 0.8)
plot(dataTa[:,0],dataTa[:,2],ls='solid', color = "red",linewidth=1,label="TC2_s_PATO", alpha = 0.8)
plot(dataTa[:,0],dataTa[:,3],ls='solid', color = "blue",linewidth=1,label="TC3_s_PATO", alpha = 0.8)
#plot(dataTa[:,0],dataTa[:,4],ls='solid', color = "yellow",linewidth=1,label="TCbottom_s", alpha = 0.8)
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
outputFile_ = "Temperature_porousMat.pdf"
savefig(outputFile_)

