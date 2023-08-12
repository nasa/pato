
########  CALCULATE PRESSURE DIFFERENCE  ###########################################################################################

mu = 1.8e-5           # dynamic viscosity
L = 0.1               # length of the tube
Q = 1e-7 *360/5       # volumetric flow rate, 5 is the angle of the wedge in degrees
pi = 3.1415926           
R = 0.005             # radius 5mm
P=8*mu*L*Q/(pi*R**4)  #Pressure difference, Hagen–Poiseuille equation

########  CALCULATE PRESSURE DIFFERENCE(SIMULATION) ###########################################################################################
from patoPlot import *
PATOFileList_=[]
PATOFileList_.append("../postProcessing/sampleDict/tube/2/line1_p.xy")          # List of data files
patoPlot = patoPlot(PATOFileList_)
DataTBas=patoPlot.dataList[0]
P_left = DataTBas[0,1]
P_right = DataTBas[1000,1]
P1 = P_left-P_right                 #Pressure difference, simulation

########  CALCULATE  ERROR BETWEEN TWO PRESSURE DIFFERENCE ###########################################################################################
P_error = abs(P-P1)/P

#print(P, P1, P_error)


########  CALCULATE VELOCITY  ###########################################################################################
G = P/L 
Umax = G*R**2/(4*mu)    #max velocity   Hagen–Poiseuille equation

#U = G/(4*mu)*(R**2-r**2)   parabolic velocity profile

########  CALCULATE VELOCITY (SIMULATION)  ###########################################################################################
from patoPlot import *
PATOFileList_=[]
PATOFileList_.append("../postProcessing/sampleDict1/tube/2/line1_U.xy")          # List of data files
patoPlot = patoPlot(PATOFileList_)
DataTBas=patoPlot.dataList[0]
U1max = DataTBas[0,1]

########  CALCULATE  ERROR BETWEEN TWO MAX VELOCITY ###########################################################################################
U_error = abs(Umax-U1max)/Umax


print (P, P1, Umax, U1max, P_error,  U_error) 

