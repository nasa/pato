# ----------------------------------------------------------------------------------------
# Loading the libraries
# ----------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import os                                 # to move through directories


# ----------------------------------------------------------------------------------------
# Data of the problem
# ----------------------------------------------------------------------------------------

mu = 1.5e-5

A = np.zeros((9,9))
b = np.zeros((9,1))

x0 = np.zeros((9,1))
x1 = np.zeros((9,1))
x2 = np.zeros((9,1))
x3 = np.zeros((9,1))
x4 = np.zeros((9,1))
x5 = np.zeros((9,1))
x6 = np.zeros((9,1))
x7 = np.zeros((9,1))
x8 = np.zeros((9,1))
x9 = np.zeros((9,1))

R = 0.24e-6
l = 1e-6
Epsilon = 0.9;

Dir = os.getcwd()              # save directory position

# ----------------------------------------------------------------------------------------
# Get results from the Ux simulation
# ----------------------------------------------------------------------------------------

os.chdir(Dir + "/Velocity_X_Direction/postProcessing/pressureDifferenceSurface/250/")      # get DeltaP in x direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_1_X = float(info_line[1])

os.chdir(Dir + "/Velocity_X_Direction/postProcessing/pressureDifferencePatch1/250/")      # get DeltaP in y direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_1_Y = float(info_line[1])

os.chdir(Dir + "/Velocity_X_Direction/postProcessing/pressureDifferencePatch2/250/")      # get DeltaP in z direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_1_Z = float(info_line[1])

os.chdir(Dir + "/Velocity_X_Direction/postProcessing/volFieldValue1/0/")      # get components of U vector

for line in open('volFieldValue.dat'):
    line = line.replace("(", " ")
    line = line.replace(")", " ")
    info_line = line.split()
    if (info_line[0]=="250"):
        U_1_X = float(info_line[1])*Epsilon
        U_1_Y = float(info_line[2])*Epsilon
        U_1_Z = float(info_line[3])*Epsilon

# ----------------------------------------------------------------------------------------
# Get results from the Uy simulation
# ----------------------------------------------------------------------------------------

os.chdir(Dir + "/Velocity_Y_Direction/postProcessing/pressureDifferenceSurface/250/")      # get DeltaP in y direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_2_Y = float(info_line[1])

os.chdir(Dir + "/Velocity_Y_Direction/postProcessing/pressureDifferencePatch1/250/")      # get DeltaP in x direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_2_X = float(info_line[1])

os.chdir(Dir + "/Velocity_Y_Direction/postProcessing/pressureDifferencePatch2/250/")      # get DeltaP in z direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_2_Z = float(info_line[1])

os.chdir(Dir + "/Velocity_Y_Direction/postProcessing/volFieldValue1/0/")      # get components of U vector

for line in open('volFieldValue.dat'):
    line = line.replace("(", " ")
    line = line.replace(")", " ")
    info_line = line.split()
    if (info_line[0]=="250"):
        U_2_X = float(info_line[1])*Epsilon
        U_2_Y = float(info_line[2])*Epsilon
        U_2_Z = float(info_line[3])*Epsilon

# ----------------------------------------------------------------------------------------
# Get results from the Uz simulation
# ----------------------------------------------------------------------------------------

os.chdir(Dir + "/Velocity_Z_Direction/postProcessing/pressureDifferenceSurface/250/")      # get DeltaP in z direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_3_Z = float(info_line[1])

os.chdir(Dir + "/Velocity_Z_Direction/postProcessing/pressureDifferencePatch1/250/")      # get DeltaP in x direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_3_X = float(info_line[1])

os.chdir(Dir + "/Velocity_Z_Direction/postProcessing/pressureDifferencePatch2/250/")      # get DeltaP in y direction

for line in open('fieldValueDelta.dat'):
    info_line = line.split()
    if (info_line[0]=="250"):
        DELTAP_3_Y = float(info_line[1])

os.chdir(Dir + "/Velocity_Z_Direction/postProcessing/volFieldValue1/0/")      # get components of U vector

for line in open('volFieldValue.dat'):
    line = line.replace("(", " ")
    line = line.replace(")", " ")
    info_line = line.split()
    if (info_line[0]=="250"):
        U_3_X = float(info_line[1])*Epsilon
        U_3_Y = float(info_line[2])*Epsilon
        U_3_Z = float(info_line[3])*Epsilon


# ----------------------------------------------------------------------------------------
# Create velocity vector times viscosity
# ----------------------------------------------------------------------------------------

b[0] = - U_1_X * mu
b[1] = - U_1_Y * mu
b[2] = - U_1_Z * mu
b[3] = - U_2_X * mu
b[4] = - U_2_Y * mu
b[5] = - U_2_Z * mu
b[6] = - U_3_X * mu
b[7] = - U_3_Y * mu
b[8] = - U_3_Z * mu


# ----------------------------------------------------------------------------------------
# Create pressure tensor
# ----------------------------------------------------------------------------------------

A = [[DELTAP_1_X / l, DELTAP_1_Y / l, DELTAP_1_Z / l, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000],
     [0.000000000000, 0.000000000000, 0.000000000000, DELTAP_1_X / l, DELTAP_1_Y / l, DELTAP_1_Z / l, 0.000000000000, 0.000000000000, 0.000000000000],  
     [0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, DELTAP_1_X / l, DELTAP_1_Y / l, DELTAP_1_Z / l], 
     [DELTAP_2_X / l, DELTAP_2_Y / l, DELTAP_2_Z / l, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000], 
     [0.000000000000, 0.000000000000, 0.000000000000, DELTAP_2_X / l, DELTAP_2_Y / l, DELTAP_2_Z / l, 0.000000000000, 0.000000000000, 0.000000000000], 
     [0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, DELTAP_2_X / l, DELTAP_2_Y / l, DELTAP_2_Z / l],  
     [DELTAP_3_X / l, DELTAP_3_Y / l, DELTAP_3_Z / l, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000],  
     [0.000000000000, 0.000000000000, 0.000000000000, DELTAP_3_X / l, DELTAP_3_Y / l, DELTAP_3_Z / l, 0.000000000000, 0.000000000000, 0.000000000000],   
     [0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, DELTAP_1_X / l, DELTAP_1_Y / l, DELTAP_3_Z / l]]


# ----------------------------------------------------------------------------------------
# Compute Permeability Tensor 
# ----------------------------------------------------------------------------------------

k = np.linalg.solve(A, b)      # evaluated as a vector: Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz
k = k/R/R                      # make it dimensionless

kxx = float(k[0])
kxy = float((k[1]+k[3])/2)
kxz = float((k[2]+k[6])/2)
kyy = float(k[4])
kyz = float((k[5]+k[7])/2)
kzz = float(k[8])


print("The components of the permeability tensor are: ")
print('Kxx =',kxx)
print('Kxy =',kxy)
print('Kxz =',kxz)
print('Kyy =',kyy)
print('Kyz =',kyz)
print('Kzz =',kzz)


# ----------------------------------------------------------------------------------------
# Diagonalize the Permeability Tensor 
# ----------------------------------------------------------------------------------------

K = [[kxx, kxy, kxz], 
     [kxy, kyy, kyz],
     [kxz, kyz, kzz]]

[evals, evecs] = np.linalg.eigh(K)
#print(evals)
#print(evecs)

v1 = evecs[:,0]
v2 = evecs[:,1]
v3 = evecs[:,2]
#print(v1 @ v2)     # to check ortogonality
#print(v1 @ v3)     # should be equal to 0
#print(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])  # should be equal to 1 

#print(np.arccos(v1[0]), np.arccos(v1[1]), np.arccos(v1[2])) # rotation angle (rad) for the first eigenvector
#print(np.arccos(v2[0]), np.arccos(v2[1]), np.arccos(v2[2])) # rotation angle (rad) for the second eigenvector
#print(np.arccos(v3[0]), np.arccos(v3[1]), np.arccos(v3[2])) # rotation angle (rad) for the third eigenvector

Diagonal = np.linalg.inv(evecs) @ K @ evecs
print("\nThe Diagonal permeability tensor is: ")
print(Diagonal)     # ovvero matrice con gli autovalori sulla diagonale

#############################################################################################################################################################################################


