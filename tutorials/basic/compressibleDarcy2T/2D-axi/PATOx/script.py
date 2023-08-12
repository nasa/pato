#!/usr/bin/env python3
# coding: utf-8

# ===========================================
# OpenFoam  (HeatTransfer Problem)
# ===========================================


# Loading the libraries
# ---------------------------------

import numpy as np
import os                                              # to move through directories
#import vtk                                             # to deal with .vtu file
#from vtk.util.numpy_support import vtk_to_numpy        # to convert values from vtk to numpy
import matplotlib.pyplot as plt                        # to plot


# Input the time folder to load
# -------------------------------------

openF_x = np.zeros(16)
openF_Tg = np.zeros(16)
openF_Ts = np.zeros(16)

Dir = os.getcwd()

index = 0
for time in range(0,16,1):
      os.chdir(Dir + "/postProcessing/singleGraph/porousMat/" + str(time) +"/")
      openF_x[index]   = np.loadtxt(open('line_p_Tg_Ta.xy','r'))[121, 0]
      openF_Tg[index] = np.loadtxt(open('line_p_Tg_Ta.xy','r'))[121, 2]
      openF_Ts[index]  = np.loadtxt(open('line_p_Tg_Ta.xy','r'))[121, 3]
      index = index + 1
#print(openF_Tg)

with open("test.xls",'wt')as f:
    for i in openF_Tg:
     print(i)


# Loading the data from .xy and .vtu
# -------------------------------------


# OpenFoam data

    
