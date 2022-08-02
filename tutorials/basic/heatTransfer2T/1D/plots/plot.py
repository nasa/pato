from patoPlot import *

# List of data files
PATOFileList_=[]
PATOFileList_.append("../postProcessing/singleGraph/15/line_T_g_T_s.xy")
PATOFileList_.append("../plots/AnalyticalT.txt")

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
ax.set_ylabel('Temperature (K)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# x limit
dataTNum=patoPlot.dataList[0]
dataTAna=patoPlot.dataList[1]
xlim = dataTNum[:,1]
# plot data
plot(dataTNum[:,0],dataTNum[:,1],"o", color = "blue",linewidth=2, alpha = 0.9)
plot(dataTNum[:,0],dataTAna[:,0],ls='solid', color = "red",linewidth=2, alpha = 0.9)
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
outputFile_ = "Tg.pdf"
savefig(outputFile_)
