# ----------------------------------------------------------------------------------------
# Loading the libraries
# ----------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import os                                 # to move through directories


# ----------------------------------------------------------------------------------------
# Data of the problem
# ----------------------------------------------------------------------------------------

mu = 2.117e-5                 # flow viscosity
R = 2.5231e-6                 # sphere radius in meter
epsilon = 0.950               # porosity 

Dir = os.getcwd()           # save directory position

# ----------------------------------------------------------------------------------------
# Get results from the Simulation at the final time = 2000
# ----------------------------------------------------------------------------------------

os.chdir(Dir + "/postProcessing/pressureDifferenceSurface/2000/")      # get DeltaP in x direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="2000"):
        DELTAP_X = float(info_line[1])

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

os.chdir(Dir + "/postProcessing/volFieldValue1/0/")      # get components of U vector

for line in open('volFieldValue.dat'):
    line = line.replace("(", " ")
    line = line.replace(")", " ")
    info_line = line.split()
    if (info_line[0]=="2000"):
        U_X = float(info_line[1])*epsilon
#        U_Y = float(info_line[2])
#        U_Z = float(info_line[3])

# ----------------------------------------------------------------------------------------
# Compute Permeability Tensor 
# ----------------------------------------------------------------------------------------

k = - mu * U_X / (DELTAP_X/ 0.00004)       # division by 40e-06 meter (= length of the domain) to compute the gradient of the pressure
k = k/R/R                                  # make it dimensionless

print("The average velocity is: ")
print("Ux =", U_X, "m/s")
print("The pressure gradient is: ")
print("gradP =", DELTAP_X/0.00004,"m/s²")  #gradient of pressure : DELTAP_X is expressed in m²/s² (Pressure divided by volumetric mass density)

print("The non-dimensional permeability is: ")
print('k =',k)

#############################################################################################################################################################################################


