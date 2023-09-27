#!/usr/bin/env python3

import patoPlot as pp

# List of data files
PATOFileList_ = []

PATOFileList_.append("../postProcessing/sampleDict/tube/2/line1_p.xy")


# Read data files
patoPlot = pp.patoPlot(PATOFileList_)

# Set up plot parameters
params = {
    'axes.labelsize': 14,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'text.usetex': False,
    'figure.figsize': [10, 5]
}
pp.rcParams.update(params)

# Temperature
# create figure
fig = pp.plt.figure()
ax = fig.add_subplot(1, 1, 1)
# set labels
ax.set_xlabel('Position (m)', fontsize=16)
ax.set_ylabel('Pressure (pa)', fontsize=16)
pp.tick_params(axis='x', top='off')
pp.tick_params(axis='y', right='off')
# x limit
DataTBas = patoPlot.dataList[0]


xlim = DataTBas[:, 1]
# plot data
pp.plot(DataTBas[:, 0], DataTBas[:, 1], ls='solid', color="blue", linewidth=2,
        alpha=0.9, label='p, Tube')


# tight layout
pp.plt.tight_layout()
# add grid
pp.plt.grid()
# put legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles, labels, loc=1, fontsize=13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "p.PNG"
pp.savefig(outputFile_)
