from patoPlot import *

# List of data files
region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../output/"+region+"/scalar/cp_plot")

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

### Heat capacity ###
# create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels
ax.set_xlabel('Temperature (K)', fontsize=16)
ax.set_ylabel('Char heat capacity (kJ/kg/K)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x and y limits
xlim(200,1800)
ylim(0.5,2.2)
# plot data
dataT=patoPlot.dataList[0]
dataCp=patoPlot.dataList[1]
plot(dataT[:,1],dataCp[:,1]*1e-3,ls='solid', color = "black",linewidth=2, alpha = 0.9)
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
outputFile_ = "cp.pdf"
savefig(outputFile_)
