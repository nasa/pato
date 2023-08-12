from patoPlot import *

# List of data files
PATOFileList_=[]
PATOFileList_.append("../postProcessing/singleGraph/cylinder/5000/line_mag(sigma)_mag(analyticalSigma).xy")
# Read data files
patoPlot = patoPlot(PATOFileList_)

# Set up plot parameters
params = {
    'axes.labelsize': 14,
    'legend.fontsize': 13,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'text.usetex': False,
    'figure.figsize': [10, 5]
}
rcParams.update(params)

# Create figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Create labels
ax.set_xlabel('Distance (m)', fontsize=16)
ax.set_ylabel('Stress magnitude (Pa)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
xlim(0.025,0.195)
ylim(7e7,2.5e8)
# plot the data
data=patoPlot.dataList[0]
plot(data[:,0],data[:,1],marker="x", markersize="10", linestyle="None", label="PATOx", alpha = 0.8)
plot(data[:,0],data[:,2],ls='solid', color = "black",linewidth=1,label="analytical", alpha = 0.8)
# tight layout
plt.tight_layout()
# add grid
plt.grid()
# put legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles,labels, loc=1)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "Stress.pdf"
savefig(outputFile_)
