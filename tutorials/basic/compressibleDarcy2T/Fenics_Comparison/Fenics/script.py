#!/usr/bin/env python3
# coding: utf-8

# ===========================================
# Heat Transfert with Advection & Convection
# ===========================================


# Loading the libraries
# ---------------------------------

from fenics import *
import numpy
import os
import shutil


# Mesh creation
# ---------------------

mesh = IntervalMesh(199,0.,0.01)


# Choice of the space functions 
# -------------------------------

P1  = FiniteElement("CG", mesh.ufl_cell(), 1)
P2  = FiniteElement("CG", mesh.ufl_cell(), 2)
P3  = FiniteElement("CG", mesh.ufl_cell(), 3)
F   = FunctionSpace(mesh, MixedElement([P1, P1, P3, P1, P1, P1]))

# Variational formulation
# -------------------------------------------------------


T_g_n_plus1, T_s_n_plus1, P_n_plus1, v_n_plus1, q_n_plus1, rho_n_plus1  = TrialFunction(F)                   # Test function
v_g, v_s, v_p, v_v, v_q, v_r = TestFunction(F)                                                               # Test function
sol_n = Function(F)   


# Definition of the boundary and initial conditions
# --------------------------------------------------------

# Initial conditions
ic_0 = Expression(('293', '323', '1e5', '0', '0', '1.1494'), degree=1)
sol_n = project(ic_0, F)


# Boundary conditions
bc_left_T  = Expression(('293'), degree=1)                            # Fixed inlet gas temperature
bc_left_P  = Expression(('100084.375'), degree=1)                     # Fixed inlet  pressure
bc_right_P = Expression(('1e5'), degree=1)                            # Fixed outlet pressure 
bc1 = DirichletBC(F.sub(0), bc_left_T,  'near(x[0],0.)')
bc2 = DirichletBC(F.sub(2), bc_right_P, 'near(x[0],0.01)')
bc3 = DirichletBC(F.sub(2), bc_left_P,  'near(x[0],0.)')

# Temporal parameters
Tend  = 0.05                                                         # End Time (s)
dt    = 1e-5                                                         # Time discretization

# Physical parameters
cp_g  = Constant(1676.)                                              # Gas heat capacity (m^2/s^2/K)
cp_s  = Constant(1000.)                                              # Solid heat capacity (m^2/s^2/K)
rho_s = Constant(300.)                                               # Solid density (kg/m^3)
h     = Constant(6000.)                                              # Heat exchange coefficient gas<-->solid (Kg/m/s^3/K)
k_s   = Constant(1.0)												 # Thermal conductivity (isotropic) (kg.m/s^3/K)

eps   = Constant(0.8)                                                # Porosity 
K     = Constant(1.6e-11)										     # Permeability in basis (x,y,z) of the mesh (m^2) 

R     = Constant(8.314471469)										 # Perfect gas constant (Kg*m^2/s^2/K/mol)
M     = Constant(28e-3)			          							 # Gas molar mass (Kg/mol)
mu    = Constant(1.35e-5)										     # Gas viscosity (kg/m/s)

# Variational formulation 

A1 = eps*sol_n[5]*cp_g/Constant(dt)*T_g_n_plus1*v_g*dx + h*T_g_n_plus1*v_g*dx + eps*cp_g*(T_g_n_plus1*sol_n[4]).dx(0)*v_g*dx                                         # Gas temperature equation
L1 = eps*sol_n[5]*cp_g/Constant(dt)*sol_n[0]*v_g*dx + h*sol_n[1]*v_g*dx         		                                                
    
A2 = (1-eps)*rho_s*cp_s/Constant(dt)*T_s_n_plus1*v_s*dx + h*T_s_n_plus1*v_s*dx + k_s*inner(grad(T_s_n_plus1),grad(v_s))*dx                                           # Solid temperature equation				   
L2 = (1-eps)*rho_s*cp_s/Constant(dt)*sol_n[1]*v_s*dx + h*sol_n[0]*v_s*dx      		        

A3 = eps*M/(R*Constant(dt))*(P_n_plus1/sol_n[0] - sol_n[2]/sol_n[0]**2*T_g_n_plus1)*v_p*dx + K*sol_n[2]*M/(mu*R*sol_n[0])*inner(grad(P_n_plus1),grad(v_p))*dx        # Pressure equation
L3 = eps*M/(R*Constant(dt))*(sol_n[2]/sol_n[0] - sol_n[2]/sol_n[0])*v_p*dx

A4 = inner(v_n_plus1,v_v)*dx                                # Velocity evaluation 
L4 = -K/(mu*eps)*sol_n[2].dx(0)*v_v*dx  

A5 = inner(q_n_plus1,v_q)*dx                                # Flow rate evaluation
L5 = sol_n[3]*sol_n[5]*v_q*dx

A6 = rho_n_plus1*v_r*dx                                     # Density evaluation
L6 = sol_n[2]*M/(R*sol_n[0])*v_r*dx

A = A1+A2+A3+A4+A5+A6
L = L1+L2+L3+L4+L5+L6

# Create the VTK files for visualization output
if (os.path.isdir('Fenics/output')):
    shutil.rmtree('Fenics/output')                                     # to delete the Output folder
os.mkdir('Fenics/output')                                          # to create the Output folder
vtkfile_T_g = File('Fenics/output/T_g.pvd')
vtkfile_T_s = File('Fenics/output/T_s.pvd')
vtkfile_P   = File('Fenics/output/P.pvd')
vtkfile_v   = File('Fenics/output/v.pvd')
vtkfile_q   = File('Fenics/output/q.pvd')
vtkfile_r   = File('Fenics/output/r.pvd')


# Temporal resolution
# ======================

t   = 0                               			 			    # Actual time

_T_g, _T_s, _P, _v, _q, _r = sol_n.split()                      # Writing the initial conditions

vtkfile_T_g << (_T_g, t)                                        
vtkfile_T_s << (_T_s, t)
vtkfile_P   << (_P, t)
vtkfile_v   << (_v, t)
vtkfile_q   << (_q, t)
vtkfile_r   << (_r, t)

ii = 0                                			   			    # Iteration
while t <= Tend :                                               # Temporal cycle
    print('temps =', t)
    
    solve(A == L, sol_n, [bc1, bc2, bc3]) 		   			    # Solving the system
    _T_g, _T_s, _P, _v, _q, _r = sol_n.split()   
    
    t  = t  + dt
    ii = ii + 1
    
    if ii%5e3 == 0:                                             # Save solution to file (VTK)
        vtkfile_T_g << (_T_g, t)
        vtkfile_T_s << (_T_s, t)
        vtkfile_P   << (_P, t)
        vtkfile_v   << (_v, t)
        vtkfile_q   << (_q, t)
        vtkfile_r   << (_r, t)

print('end of the temporal cycle')        
