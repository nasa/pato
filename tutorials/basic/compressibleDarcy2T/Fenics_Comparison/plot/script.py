#!/usr/bin/env python3
# coding: utf-8

# ===========================================
# OpenFoam vs Fenics (HeatTransfer Problem)
# ===========================================


# Loading the libraries
# ---------------------------------

import numpy as np
import os                                              # to move through directories
import vtk                                             # to deal with .vtu file
from vtk.util.numpy_support import vtk_to_numpy        # to convert values from vtk to numpy
import matplotlib.pyplot as plt                        # to plot


# Input the time folder to load
# -------------------------------------

time = 0.05


# Loading the data from .xy and .vtu
# -------------------------------------

Dir = os.getcwd()

fenics_file = int(time/0.05)
if (fenics_file<10):
    fenics_file = "00000" + str(fenics_file)
elif (fenics_file<100):
    fenics_file = "0000" + str(fenics_file)
elif (fenics_file<1000):
    fenics_file = "000" + str(fenics_file)
elif (fenics_file<10000):
    fenics_file = "00" + str(fenics_file)
elif (fenics_file<100000):
    fenics_file = "0" + str(fenics_file)
elif (fenics_file<1000000):
    fenics_file = str(fenics_file)

# OpenFoam data

os.chdir(Dir + "/postProcessing/singleGraph/" + str(time) + "/")
if (time>0):
    openF_x   = np.loadtxt(open('line_T_g_T_s_p_rho_g.xy','r'))[:, 0]
    openF_x = np.concatenate(([0], openF_x))
    openF_Tg  = np.loadtxt(open('line_T_g_T_s_p_rho_g.xy','r'))[:, 1]
    openF_Tg = np.concatenate(([293], openF_Tg))
    openF_Ts  = np.loadtxt(open('line_T_g_T_s_p_rho_g.xy','r'))[:, 2]
    openF_Ts = np.concatenate(([openF_Ts[0]], openF_Ts))
    openF_P   = np.loadtxt(open('line_T_g_T_s_p_rho_g.xy','r'))[:, 3]
    openF_P = np.concatenate(([100084.375], openF_P))
    openF_rho = np.loadtxt(open('line_T_g_T_s_p_rho_g.xy','r'))[:, 4]
    openF_rho = np.concatenate(([1.1494], openF_rho))
else:
    openF_x   = np.loadtxt(open('line_T_g_T_s_p.xy','r'))[:, 0]
    openF_x   = np.concatenate(([0], openF_x))
    openF_Tg  = np.loadtxt(open('line_T_g_T_s_p.xy','r'))[:, 1]
    openF_Tg  = np.concatenate(([293], openF_Tg))
    openF_Ts  = np.loadtxt(open('line_T_g_T_s_p.xy','r'))[:, 2]
    openF_Ts  = np.concatenate(([openF_Ts[0]], openF_Ts))
    openF_P   = np.loadtxt(open('line_T_g_T_s_p.xy','r'))[:, 3]
    openF_P   = np.concatenate(([openF_P[0]], openF_P))
    openF_rho = np.zeros(len(openF_P))
    
openF_q   = np.loadtxt(open('line_u_g_q.xy','r'))[:, 4]
openF_q = np.concatenate(([openF_q[0]], openF_q))
openF_v   = np.loadtxt(open('line_u_g_q.xy','r'))[:, 1]
openF_v = np.concatenate(([openF_v[0]], openF_v))

# Fenics data

reader = vtk.vtkXMLUnstructuredGridReader()
os.chdir(Dir + "/Fenics/output/")
reader.SetFileName('T_g'+fenics_file+'.vtu')
reader.Update()  
fenics_x  = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())[:,0]             
fenics_Tg = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars())

reader.SetFileName('T_s'+fenics_file+'.vtu')
reader.Update()  
fenics_Ts = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars())

reader.SetFileName('P'+fenics_file+'.vtu')
reader.Update()  
fenics_P = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars())

reader.SetFileName('r'+fenics_file+'.vtu')
reader.Update()  
fenics_rho = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars())

reader.SetFileName('q'+fenics_file+'.vtu')
reader.Update()  
fenics_q = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars())

reader.SetFileName('v'+fenics_file+'.vtu')
reader.Update()  
fenics_v = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars())


# plot the comparisons
# -------------------------------------

if not os.path.isdir(Dir + "/plot/figures"):
    os.mkdir(Dir + "/plot/figures")
os.chdir(Dir + "/plot/figures")

# Temperature

plt.plot(openF_x,openF_Tg,'b', linewidth=1, label="GAS")
plt.plot(openF_x,openF_Ts,'k', linewidth=1, label="SOLID")
plt.plot(fenics_x[0:-1:7],fenics_Tg[0:-1:7],'bo', markersize=7, mfc='none') 
plt.plot(fenics_x[0:-1:7],fenics_Ts[0:-1:7],'ko', markersize=7, mfc='none', label="fenics")
plt.plot(openF_x,openF_Ts,'k', linewidth=1, label="OpenFoam")
plt.title("Time t = " + str(time))
plt.xlabel ("x [m]");
plt.ylabel ("Temperature [k]");
plt.legend(loc='lower right', borderaxespad=0., frameon = False)
plt.savefig('temperatures_'+str(time)+'.pdf')
plt.show()

# Pressure

plt.plot(openF_x,openF_P,'k', linewidth=1, label="OpenFoam")
plt.plot(fenics_x[0:-1:7],fenics_P[0:-1:7],'ko', markersize=7, mfc='none', label="fenics") 
plt.title("Time t = " + str(time))
plt.xlabel ("x [m]");
plt.ylabel ("Pressure [Pa]");
plt.legend(loc='lower right', borderaxespad=0., frameon = False)
plt.savefig('pressure_'+str(time)+'.pdf')
plt.show()

# Density

plt.plot(openF_x,openF_rho,'k', linewidth=1, label="OpenFoam")
plt.plot(fenics_x[0:-1:7],fenics_rho[0:-1:7],'ko', markersize=7, mfc='none', label="fenics") 
plt.title("Time t = " + str(time))
plt.xlabel ("x [m]");
plt.ylabel ("Density [Kg/m^3]");
plt.legend(loc='lower right', borderaxespad=0., frameon = False)
plt.savefig('density_'+str(time)+'.pdf')
plt.show()

# Flow rate

plt.plot(openF_x,openF_q,'k', linewidth=1, label="OpenFoam")
plt.plot(fenics_x[0:-1:7],fenics_q[0:-1:7],'ko', markersize=7, mfc='none', label="fenics") 
plt.title("Time t = " + str(time))
plt.xlabel ("x [m]");
plt.ylabel ("Flow rate [Kg/m^2/s]");
plt.legend(loc='lower right', borderaxespad=0., frameon = False)
plt.savefig('flow_'+str(time)+'.pdf')
plt.show()

# Velocity

plt.plot(openF_x,openF_v,'k', linewidth=1, label="OpenFoam")
plt.plot(fenics_x[0:-1:7],fenics_v[0:-1:7],'ko', markersize=7, mfc='none', label="fenics") 
plt.title("Time t = " + str(time))
plt.xlabel ("x [m]");
plt.ylabel ("Velocity [m/s]");
plt.legend(loc='lower right', borderaxespad=0., frameon = False)
plt.savefig('velocity_'+str(time)+'.pdf')
plt.show()

