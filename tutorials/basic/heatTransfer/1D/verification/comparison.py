#########################################################################################################################

# python code to compare the analytical and finite difference solutions
# with the openFoam ones 

#########################################################################################################################

########    LIBRARIES    ################################################################################################

import numpy as np
import os                                              # to move through directories
import matplotlib.pyplot as plt                        # to plot

########    TIME TO COMPARE   ###########################################################################################

time = 60                                                                       

########    OPEN FOAM LOADING   #########################################################################################

Dir = os.getcwd()                                                               # current directory

os.chdir(Dir + "/verification/postProcessing/singleGraph/" + str(time) + "/")

openFoam_x = np.loadtxt(open('line_T.xy','r'))[:, 0]
openFoam_T = np.loadtxt(open('line_T.xy','r'))[:, 1]

########    FINITE DIFFERENCE LOADING   #################################################################################

os.chdir(Dir + "/verification/analytical_and_finiteDifference/")

with open("finite_difference.txt") as f:
   i = 0
   for line in f:                                             # to count the lines of the file
      i = i + 1 

finiteDifference_x = np.zeros(i)
finiteDifference_T = np.zeros(i)

i = 0
with open("finite_difference.txt") as f:
    for line in f:
      s = line.split()
      finiteDifference_x[i] = s[0] 
      finiteDifference_T[i] = s[int(time/1)] 
      i = i + 1

########    ANALYTICAL SOLUTION LOADING   ################################################################################

os.chdir(Dir + "/verification/analytical_and_finiteDifference/")

with open("analytical_solution.txt") as f:
   i = 0
   for line in f:                                             # to count the lines of the file
      i = i + 1 

analyticalSolution_x = np.zeros(i)
analyticalSolution_T = np.zeros(i)

i = 0
with open("analytical_solution.txt") as f:
    for line in f:
      s = line.split()
      analyticalSolution_x[i] = s[0] 
      analyticalSolution_T[i] = s[int(time/1)] 
      i = i + 1

########    PLOT   ######################################################################################################

plt.figure()
plt.plot(analyticalSolution_x,analyticalSolution_T, 'k', linewidth=0.7, label="Analytical Solution")
plt.plot(openFoam_x[0:-1:2],openFoam_T[0:-1:2], 'ro', markersize=4, mfc='none', label="OpenFOAM Solution")
plt.plot(finiteDifference_x[0:-1:2],finiteDifference_T[0:-1:2], 'bo', markersize=4, mfc='none', label="Finite Difference Solution")
plt.title('1D Heat Problem')
plt.xlabel('X position [m]')
plt.ylabel('Temperature [k]')
plt.legend(loc='upper left', borderaxespad=1., frameon = False)
plt.savefig('Comparison.pdf')
plt.show()




