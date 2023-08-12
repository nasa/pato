from patoPlot import *

# List of data files
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../output/"+region+"/scalar/RR.O2_plot")
PATOFileList_.append("../output/"+region+"/scalar/RR.CO_plot")
PATOFileList_.append("../output/"+region+"/scalar/RR.C(gr)_plot")
PATOFileList_.append("../output/"+region+"/scalar/RR.CO2_plot")

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

### Reaction rates ###
# create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels
ax.set_xlabel('Temperature (K)', fontsize=16)
ax.set_ylabel('Reaction rate per specie (mg/m3/s)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limits
xlim(300,3000)
# plot data
dataT=patoPlot.dataList[0]
dataRR=patoPlot.dataList[1]
plot(dataT[:,1],dataRR[:,1]*1e6,ls='solid', color = "red",linewidth=2, label= "O2", alpha = 0.9)
dataRR=patoPlot.dataList[3]
plot(dataT[:,1],dataRR[:,1]*1e6,ls='dashed', color = "red",linewidth=2, label = "C(gr)", alpha = 0.9)
dataRR=patoPlot.dataList[2]
plot(dataT[:,1],dataRR[:,1]*1e6,ls='dashed', color = "green",linewidth=2, label = "CO", alpha = 0.9)
dataRR=patoPlot.dataList[4]
plot(dataT[:,1],dataRR[:,1]*1e6,ls='solid', color = "green",linewidth=2, label = "CO2", alpha = 0.9)
# tight layout
plt.tight_layout()
# add grid
plt.grid()
# put legends
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=2, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "RR.pdf"
savefig(outputFile_)
