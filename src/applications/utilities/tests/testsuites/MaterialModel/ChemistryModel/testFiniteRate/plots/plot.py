from pylab import *
import brewer2mpl
from numpy import *
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')

region = "porousMat"
PATOFileList_=[]
PATOFileList_.append("../output/"+region+"/scalar/Ta_plot")
PATOFileList_.append("../output/"+region+"/scalar/h_g_plot")

# brewer2mpl.get_map args: set name  set type  number of colors
bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors

params = {
    'axes.labelsize': 14,
    'text.fontsize': 14,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'text.usetex': False,
    'figure.figsize': [10, 5]
}
rcParams.update(params)

def num(s):
    table=[]
    try:
        for i in s:
            table.append(float(i))
        return table
    except ValueError:
        return table

def readDataPATO(fileName):
    data = []
    title = []
    with open(fileName) as f:
        content = f.readlines()
        i = 0
        for line in content:
            if i == 0:
                title=line.split()
            if i > 0:
                values = line.split()
                data.append(num(values))
            i = i + 1
    data=filter(None,data)
    first = len(data[0])
    addColumn = 1
    for row in data:
        lenRow = len(row)
        if lenRow<first:
            for i in range(0,first-lenRow):
                row.insert(i+addColumn,None)

    data=array(data)
    return data,title

dataList=[]
titleList = []
fileNameList = PATOFileList_
lsList = ["solid"]

for fileName in fileNameList:
    data,title= readDataPATO(fileName)
    dataList.append(data)
    titleList.append(title)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('Temperature (K)', fontsize=16)
ax.set_ylabel('Gas enthalpy (MJ/kg)', fontsize=16)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
xlim(300,1800)
dataT=dataList[0]
dataH=dataList[1]
plot(dataT[:,1],dataH[:,1]/1e6,ls='solid', color = "black",linewidth=2, alpha = 0.9)
handles, labels = ax.get_legend_handles_labels()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')
outputFile_ = "h_g.pdf"
savefig(outputFile_)
