from patoPlot import *

species=[
        "H2"
       ,"CH4"
       ,"CO"
       ,"CO2"
       ,"A1OH"
       ,"H2O"]

speciesName=["H2"
        ,"CH4"
       ,"CO"
       ,"CO2"
       ,"A1OH"
       ,"H2O"]

sizeSpecies=len(species)
PATOFileList_ = []

for specieI in species:
        PATOFileList_.append("../output/porousMat/scalar/normalizedPi["+specieI+"]_plot")
for specieI in species:
    PATOFileList_.append("../data-3.1/" + specieI + ".txt")
PATOFileList_.append("../output/porousMat/scalar/normalizedPiTotal_plot")
PATOFileList_.append("../data-3.1/normalizedPiTotal.txt")
PATOFileList_.append("../output/porousMat/scalar/rhoRatio_plot")
PATOFileList_.append("../data-3.1/rhoRatio.txt")
PATOFileList_.append("../output/porousMat/scalar/Ta_plot")
PATOFileList_.append("../data-3.1/Temperature.txt")

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
lsList = ["solid"]


temperaturePATO = patoPlot_.dataList[len(PATOFileList_) - 2][:, 1]
temperatureData = patoPlot_.dataList[len(PATOFileList_) - 1][1:, 1]

rhoRatioPATO = patoPlot_.dataList[len(PATOFileList_) - 4][:, 1]
rhoRatioData = patoPlot_.dataList[len(PATOFileList_) - 3][:, 1]
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('$T (K)$', fontsize=16)
ax.set_ylabel('$\\rho_s/\\rho_{s,v}$ $(-)$', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
xlim(300, 1463.35)
ylim(0.75,1)
plot(temperaturePATO,rhoRatioPATO,ls='solid', color = "black",linewidth=2,label="PATO", alpha = 0.8)
plot(temperatureData,rhoRatioData,"ok",ms=3,label="data", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =16)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "rhoRatio.pdf"
savefig(outputFile_)

normalizedPiTotalPATO = patoPlot_.dataList[len(PATOFileList_) - 6][:, 1]
normalizedPiTotalData = patoPlot_.dataList[len(PATOFileList_) - 5][:, 1]
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('$T (K)$', fontsize=16)
ax.set_ylabel('$\pi_{tot}/\\rho_{s,v}$ $(10^{-3}/s)$', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
xlim(300, 1463.35)
plot(temperaturePATO,normalizedPiTotalPATO*1e3,ls='solid', color = "black",linewidth=2,label="PATO: total", alpha = 0.8)
data = patoPlot_.dataList[0]
plot(temperaturePATO,data[:, 1] * 1e3,"--k",linewidth=1,label="PATO: reaction 1", alpha = 0.8)
data = patoPlot_.dataList[1]
plot(temperaturePATO,data[:, 1] * 1e3,":k",linewidth=1,label="PATO: reaction 2", alpha = 0.8)
data = patoPlot_.dataList[2]
plot(temperaturePATO,data[:, 1] * 1e3,"-.k",linewidth=1,label="PATO: reaction 3", alpha = 0.8)
data = patoPlot_.dataList[3]
plot(temperaturePATO,data[:, 1] * 1e3,"-k",linewidth=1,label="PATO: reaction 4", alpha = 0.8)
plot(temperatureData,normalizedPiTotalData*1e3,"ok",ms=3,label="data", alpha = 0.8)
handles, labels = ax.get_legend_handles_labels()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =16)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "normalizedPiTotal.pdf"
savefig(outputFile_)
