from patoPlot import *

PATOFileList_=[]
PATOFileList_.append("../output/porousMat/scalar/Ta_plot")
PATOFileList_.append("../output/porousMat/scalar/Ta_surfacePatch")
PATOFileList_.append("../data/Ta_ref")

patoPlot_ = patoPlot(PATOFileList_)

params = {
    'axes.labelsize': 14,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'text.usetex': False,
    'figure.figsize': [10, 5]
}
rcParams.update(params)

fileNameList = PATOFileList_

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('$t (s)$', fontsize=16)
ax.set_ylabel('$T$ $(K)$', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
xlim(1,120)
data=patoPlot_.dataList[1]
plot(data[:,0],data[:,1],ls='solid', color = "black",linewidth=2,label="PATO-Tw", alpha = 0.8)
data=patoPlot_.dataList[0]
plot(data[:,0],data[:,1],ls='solid', color = "purple",linewidth=2,label="PATO-TC1", alpha = 0.8)
plot(data[:,0],data[:,2],ls='solid', color = "red",linewidth=2,label="PATO-TC2", alpha = 0.8)
plot(data[:,0],data[:,3],ls='solid', color = "blue",linewidth=2,label="PATO-TC3", alpha = 0.8)
plot(data[:,0],data[:,4],ls='solid', color = "green",linewidth=2,label="PATO-TC4", alpha = 0.8)
data=patoPlot_.dataList[2]
plot(data[:,0],data[:,1],ls='dashed', color = "purple",linewidth=2,label="ref-TC1", alpha = 0.8)
plot(data[:,0],data[:,2],ls='dashed', color = "red",linewidth=2,label="ref-TC2", alpha = 0.8)
plot(data[:,0],data[:,3],ls='dashed', color = "blue",linewidth=2,label="ref-TC3", alpha = 0.8)
plot(data[:,0],data[:,4],ls='dashed', color = "green",linewidth=2,label="ref-TC4", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =11)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "Ta.pdf"
savefig(outputFile_)

