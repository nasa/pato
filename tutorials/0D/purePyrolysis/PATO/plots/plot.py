from patoPlot import *

# List of data files
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../output/"+region+"/scalar/piTotal_plot")

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

### Total pyrolysis production rates ###
# create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels
ax.set_xlabel('Temperature (K)', fontsize=16)
ax.set_ylabel('Total production rates (kg/m3/s)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# set limits
xlim(300,1388)
# plot data
dataT=patoPlot.dataList[0]
dataTot=patoPlot.dataList[1]
plot(dataT[:,1],dataTot[:,1],ls='solid', color = "black",linewidth=2, alpha = 0.9)
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
outputFile_ = "piTot.pdf"
savefig(outputFile_)
