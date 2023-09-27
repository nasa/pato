# ----------------------------------------------------------------------------------------
# Loading the libraries
# ----------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import os                                 # to move through directories
from math import pi

# ----------------------------------------------------------------------------------------
# Data of the problem
# ----------------------------------------------------------------------------------------

mu = 1.656e-5                 # flow viscosity
R = 8.003e-6                 # sphere radius in meter
epsilon = 0.497              # porosity 
characteristic_length = 13.115e-6
Temperature = 300            #Kelvins
molar_weight = 28.0134e-3    #kg
nbBufferCells = 100
nbDomainCells = 200
nbTotalCells = 31744
Dir = os.getcwd()           # save directory position

# ----------------------------------------------------------------------------------------
# Get results from the Simulation at the final time = 2e-4
# ----------------------------------------------------------------------------------------

os.chdir(Dir + "/0.0002/")      # get DeltaP in x direction +remettre /Simulation au début

Lux = [] 
LV = []                        #  Initialisation
S =0
VolumeTot = 0
counter1 = 0
counter2 = 0
counter3 = 0

for line in open('UMean.txt'):
    info_line = line.split()
    if counter1 >= 22 and counter1 <=(nbTotalCells + 23): 
        if len(info_line) == 3:
            info0 = info_line[0]
            info2 = info_line[2]
            if info0[:1] == '(' :
                endL = len(info_line[2]) - 1
                if info2[endL:] == ')':
                    if info0 != '(' and info0 !=  ')' and info0 !=  '{' and info0 !=  '}' and info0 !=  ';' and info0 != 'boundaryField'and info0 != 'inlet'and info0 != 'outlet'and info0 != 'wall_buffer_1' and info0 != 'wall_buffer_2' and info0 != 'wall_buffer_3' and info0 != 'wall_buffer_4' and info0 != 'wall_domain_1' and info0 != 'wall_domain_2' and info0 != 'z1' and info0 != 'z2'and info0 != 'solid1'and info0 != 'solid2'and info0[:1] != '/':
                        ux = float(info_line[0].lstrip('('))
                        Lux = Lux + [ux]
    counter1 = counter1 + 1

print('length of UMean list = ', len(Lux))

if len(Lux) == nbTotalCells:
    print("OK !")
else:
    print("Problem with reading of file UMean !")

for Line in open('V.txt'):
    info_Line = Line.split()
    if counter2 >= 22:
        if len(info_Line) == 1:
            Info0 = info_Line[0]
            if Info0 != '(' and Info0 !=  ')' and Info0 !=  '{' and Info0 !=  '}' and Info0 !=  ';' and Info0 != 'boundaryField'and Info0 != 'inlet'and Info0 != 'outlet'and Info0 != 'wall_buffer_1' and Info0 != 'wall_buffer_2' and Info0 != 'wall_buffer_3' and Info0 != 'wall_buffer_4' and Info0 != 'wall_domain_1' and Info0 != 'wall_domain_2' and Info0 != 'z1' and Info0 != 'z2'and Info0 != 'solid1'and Info0 != 'solid2'and Info0[:1] != '/':
                v = float(Info0)
                LV = LV + [v]
    counter2 = counter2 + 1

print( 'length of volume list = ' , len(LV))

if len(LV) == nbTotalCells:
    print("OK !")
else:
    print("Problem with reading of file V !")

for u in Lux :
    if counter3 >= (nbTotalCells//2) and counter3 <= (nbTotalCells - (nbBufferCells*nbBufferCells)):
        VolCell = LV[counter3]
        uCell = u
        S = S + uCell*VolCell
        VolumeTot = VolumeTot + VolCell
    counter3 = counter3 + 1

U_X = S / VolumeTot

# Knudsen number calculation: 

counter4 = 0
counter5 = 0
RHONMean = 0

for line_rhoNMean in open('rhoNMean.txt'):
    info_line_rhoNMean = line_rhoNMean.split()
    if  counter4 >= 22 and counter4 <= 6250:
        if len (info_line_rhoNMean) == 1:
            info_rhoNMean = info_line_rhoNMean[0]
            if info_rhoNMean != '(' and info_rhoNMean !=  ')' and info_rhoNMean !=  '{' and info_rhoNMean !=  '}' and info_rhoNMean !=  ';' and info_rhoNMean != 'boundaryField'and info_rhoNMean != 'inlet'and info_rhoNMean != 'outlet'and info_rhoNMean != 'wall_buffer_1' and info_rhoNMean != 'wall_buffer_2' and info_rhoNMean != 'wall_buffer_3' and info_rhoNMean != 'wall_buffer_4' and info_rhoNMean != 'wall_domain_1' and info_rhoNMean != 'wall_domain_2' and info_rhoNMean != 'z1' and info_rhoNMean != 'z2'and info_rhoNMean != 'solid1'and info_rhoNMean != 'solid2':
                 rhoNMean = float(info_line_rhoNMean[0])
                 RHONMean = RHONMean + rhoNMean
                 counter5 = counter5 + 1
    counter4 = counter4 + 1

Ndensity = RHONMean / counter5
print ('Ndensity = ' , Ndensity)
rho = Ndensity *(molar_weight/6.02e23)
Pressure = rho *8.314 * Temperature /(molar_weight)
mean_free_path = (9.5e-8) * (1e5) *Temperature /(298 * Pressure)
Kn = mean_free_path/ characteristic_length
print ('Kn = ', Kn)

# Determination of the variation of pressure

os.chdir(Dir + "/postProcessing/pressureDifferenceSurface/0/")      # get DeltaP in x direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="0.0002"):
        DELTAP_X = float(info_line[1])

#os.chdir(Dir + "/postProcessing/pressureDifferenceSurface/0/")      # get DeltaP in x direction +remettre /Simulation au début

#for line in open('fieldValueDelta.dat'):
#    info_line = line.split()
#    if (info_line[0]=="0.00012"):
#        DELTAP_X = float(info_line[1])

#os.chdir(Dir + "/Velocity_X_Direction/postProcessing/pressureDifferencePatch1/2000/")      # get DeltaP in y direction
#
#for line in open('fieldValueDelta.dat'):
#    info_line = line.split()
#    if (info_line[0]=="2000"):
#        DELTAP_Y = float(info_line[1])
#
#os.chdir(Dir + "/Velocity_X_Direction/postProcessing/pressureDifferencePatch2/2000/")      # get DeltaP in z direction
#
#for line in open('fieldValueDelta.dat'):
#    info_line = line.split()
#    if (info_line[0]=="2000"):
#        DELTAP_Z = float(info_line[1])

#os.chdir(Dir + "/postProcessing/volFieldValue1/0/")      # get components of U vector +remettre /Simulation au début

#for line in open('volFieldValue.dat'):
#    line = line.replace("(", " ")
#    line = line.replace(")", " ")
#    info_line = line.split()
#    if (info_line[0]=="0.00012"):
#        U_X = float(info_line[1])*epsilon
#        U_Y = float(info_line[2])
#        U_Z = float(info_line[3])

# ----------------------------------------------------------------------------------------
# Compute Permeability Tensor 
# ----------------------------------------------------------------------------------------

k = - mu * U_X / (DELTAP_X/0.00002)       # division by 20e-06 meter (= length of the studied domain - around the second cylinder) to compute the gradient of the pressure
k = k/(2*R*2*R)                                  # make it dimensionless

print("The average velocity is: ")
print("Ux =", U_X, "m/s")
print("DELTA_P = ", DELTAP_X)
print("The pressure gradient is: ")
print("gradP =", DELTAP_X/0.00002,"m/s²")  #gradient of pressure : DELTAP_X is expressed in m²/s² (Pressure divided by volumetric mass density)

print("The non-dimensional permeability is: ")
print('k =',k)

#############################################################################################################################################################################################


