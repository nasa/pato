#########################################################################################################################

# python code for solving the 1D heat problem:
#   dt(T) = alpha * dxx(T)
#       with Dirichlet b.c. 
#           T(t,xl) = T0
#           T(t,xr) = Tu
#       and initial condition
#           T(0,x) = T0

#########################################################################################################################

def analytical_solution(z,t):
    a = (Tu-T0)/L * z
    summation = 0
    for n in range(1,100):
        b = 2*(Tu-T0)/(np.pi*n) * (-1)**n
        c = np.sin(n*np.pi*z/L)
        d = np.exp(-k/(rho*Cp)*n**2*np.pi**2/L**2*t)
        summation = summation + b*c*d
    return (T0 + a + summation)

########    LIBRARIES    ################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import os                                              # to move through directories

########    DATE    ######################################################################################################

L = 0.05                        # length of the boundary

Tu = 573                        # imposed temperature at t=0 
T0 = 273                        # initial temperature 

nx = 50                         # space discretization                  
dt = 1e-3                       # time  discretization

x = np.linspace(0,L,nx)     

k = 1                           # thermal conductivity 
rho = 1000                      # density
Cp  = 1000                      # specific heat at P constant
alpha = k/(rho*Cp)              # thermal diffusivity
     
time = np.linspace(1,60,60)     # saving times

########    ANALYTICAL SOLUTION    ########################################################################################

Temp_analytical = np.zeros((nx,len(time)))  

for i in range(0,nx):
    for j in range(0,len(time)):
        Temp_analytical[i,j] = analytical_solution(x[i],time[j])    

########    NUMERICAL SOLUTION    #########################################################################################

dx = x[1] - x[0]

Temp_numerical = np.zeros((nx,len(time)))

Temperature = T0*np.ones(nx)                    # initial condition
Temperature[-1] = Tu                            # boundary condition

t = 0
conta = 0

while (t<=time[-1]):      
    while (t<=time[conta]):                      
        T_old = Temperature 
        for i in range (1,nx-1):                                                 
             Temperature[i] = T_old[i] + alpha * dt / dx**2 * ( (T_old[i-1]+T_old[i+1]-2*T_old[i]) )         # finite difference solution
        t = t + dt
        print('time = ', t, '[s]')
    Temp_numerical[:,conta] = Temperature  
    conta = conta + 1

########    SAVE TO FILE    ################################################################################################

os.chdir(os.getcwd() + "/verification/analytical_and_finiteDifference/")

f = open('finite_difference.txt', 'w')

for i in range (0,nx):
    f.write('%f' %x[i])
    f.write("\t")
    for j in range (0,len(time)):
        f.write('%f' %Temp_numerical[i,j])
        f.write("\t")
    f.write("\n")

f.close()

f = open('analytical_solution.txt', 'w')

for i in range (0,nx):
    f.write('%f' %x[i])
    f.write("\t")
    for j in range (0,len(time)):
        f.write('%f' %Temp_analytical[i,j])
        f.write("\t")
    f.write("\n")

f.close()

########    PLOT    #######################################################################################################

plt.figure()
plt.clf()
plt.plot(x,Temp_analytical[:,0],  'k',  linewidth=0.7, label="Analytical Solution")
plt.plot(x,Temp_analytical[:,4], 'k',  linewidth=0.7)
plt.plot(x,Temp_analytical[:,10], 'k',  linewidth=0.7)
plt.plot(x,Temp_analytical[:,25], 'k',  linewidth=0.7)
plt.plot(x,Temp_analytical[:,40], 'k',  linewidth=0.7)
plt.plot(x,Temp_analytical[:,-1], 'k',  linewidth=0.7)
plt.plot(x,Temp_numerical[:,0],  'ko', markersize=4, mfc='none', label="Finite Difference Solution")
plt.plot(x,Temp_numerical[:,4], 'ko', markersize=4, mfc='none')
plt.plot(x,Temp_numerical[:,10], 'ko', markersize=4, mfc='none')
plt.plot(x,Temp_numerical[:,25], 'ko', markersize=4, mfc='none')
plt.plot(x,Temp_numerical[:,40], 'ko', markersize=4, mfc='none')
plt.plot(x,Temp_numerical[:,-1], 'ko', markersize=4, mfc='none')
plt.title('1D Heat Problem')
plt.xlabel('X position [m]')
plt.ylabel('Temperature [k]')
plt.legend(loc='upper left', borderaxespad=1., frameon = False)
arrow_properties = dict(facecolor="black", width=0.5, headwidth=4, shrink=0.1)
plt.annotate("", xy=(0.4*x[-1], 0.82*Tu), xytext=(1.05*x[-1], 0.5*Tu), arrowprops=arrow_properties)
plt.text(0.40*x[-1], 0.805*Tu, 'time', fontsize=11)
plt.savefig('Analytical_and_finiteDifference.pdf')
plt.show()
