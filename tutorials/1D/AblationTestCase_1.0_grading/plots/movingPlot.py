import matplotlib.pyplot as plt
import numpy as np
from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy
from pylab import *
import brewer2mpl
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')


#Read data

times = range(121)
PATOFileList_=[]
for i in times:
    PATOFileList_.append("../output/porousMat/oxidation/rho_s_" + str(i))

for i in times:
    PATOFileList_.append("../../AblationTestCase_1.0/output/porousMat/oxidation/rho_s_" + str(i))

# brewer2mpl.get_map args: set name  set type  number of colors
bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors

params = {
    'axes.labelsize': 14,
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

# fig =  plt.figure()
fig, ax = plt.subplots(1,figsize=(8,5), facecolor='white')
# ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('in-depth (cm)', fontsize=16)
ax.set_ylabel('solid density (kg/m3)', fontsize=16)
xlim(0,2)
ylim(200,400)
tick_params(axis='x', top='off')
tick_params(axis='y', right='off')
# data=dataList[0]
xx = (.05-dataList[0][:,0])*100
zz = lambda t:dataList[int(t)][:,1]
zz2 = lambda t:dataList[int(t)+len(times)][:,1]
ax.plot(xx,zz(0), ls='solid', color = "red",linewidth=2, label = "graded TACOT, t=0s", alpha = 0.9)
ax.plot(xx,zz2(0), ls='dashed', color = "red",linewidth=2, label = "TACOT, t=0s", alpha = 0.9)
line, = ax.plot(xx,zz(0), ls='solid', color = "black",linewidth=2, alpha = 0.9)
line2, = ax.plot(xx,zz2(0), ls='dashed', color = "black",linewidth=2, alpha = 0.9)
handles, labels = ax.get_legend_handles_labels()
plt.grid()
legend = ax.legend(handles,labels, loc=1, fontsize =13)
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')



# DRAW A FIGURE WITH MATPLOTLIB

duration = 6
fps = 20
# ANIMATE WITH MOVIEPY (UPDATE THE CURVE FOR EACH t). MAKE A GIF.

def make_frame_mpl(t):
    line.set_ydata( zz(fps*t))  # <= Update the curve
    line.set_label("graded TACOT, t="+str(int(fps*t))+"s")
    line2.set_ydata(zz2(fps * t))  # <= Update the curve
    line2.set_label("TACOT, t="+str(int(fps*t))+"s")
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles,labels, loc=1, fontsize =13)
    return mplfig_to_npimage(fig) # RGB image of the figure

animation =mpy.VideoClip(make_frame_mpl, duration=duration)
animation.write_gif("rho_s.gif", fps=fps)
