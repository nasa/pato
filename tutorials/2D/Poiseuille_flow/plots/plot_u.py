from patoPlot import *

# List of data files
PATOFileList_=[]

PATOFileList_.append("../postProcessing/sampleDict1/tube/2/line1_U.xy")
PATOFileList_.append("./line_u1.txt")

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
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Velocity (m/s)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
DataTBas=patoPlot.dataList[0]
DataT1=patoPlot.dataList[1]

xlim = DataTBas[:,1]
# plot data
plot(DataTBas[:,0],DataTBas[:,1],"o", color = "black",linewidth=2, alpha = 0.9,label='u_simulation, y direction')
plot(DataTBas[:,0],DataT1[:,1],ls='solid', color = "red",linewidth=2, alpha = 0.9 ,label='u_analytical, y direction')

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
outputFile_ = "u.PNG"
savefig(outputFile_)
