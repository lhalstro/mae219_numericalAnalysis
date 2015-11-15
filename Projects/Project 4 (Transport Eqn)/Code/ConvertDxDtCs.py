import numpy as np
import matplotlib.pyplot as plt 
import os
import time
import pandas as pd
from scipy import stats

global u, D, L, k, tend
u = 0.2                 #convection velocity [m/s]
D = 0.005               #diffusion coefficient [m^2/s]
L = 1.                  #length [m]
k = 2*np.pi/L           #constant [m^-1]
tend = 1./(k**2 * D)     #tau (max time) [s]

def dxdt2Cs (dx, dt):
    C, s = u * dt/dx, D * dt/dx**2
    return C, s

def Cs2dxdt(C, s):
    dx = (D*C) / (u*s)
    dt = s * dx**2 / D
    return dx, dt

def StabilityCriteria(name, C, s, dx, dt):
    """Assess stability of given numerical method for given case
    f --> array of data for given method
    dx --> grid spacing for given case
    dt --> time spacing for given case
    """
    print('C=%1.4f, s=%1.4f, dx=%1.4f, dt=%1.4f'%(C,s,dx,dt))
    #SELECT STABILITY CRITERIA FOR SPECIFIC METHOD
    if name=='FTCS' and dx<=2*D/u and dt<=dx**2/(2*D):
        #FTCS stability criterion
        print('FTCS stable')
    elif name=='Upwind' and C+2*s<1:
        #Upwind stability criterion (CFL)
            print('upwind stable')
    elif name=='Trapezoidal':
        #Trapezoidal is implicit, thus inhernetly stable
        print('Trapezoidal stable')
    elif name=='QUICK' and C<=min(2-4*s, (2*s)**0.5):
        #QUICK stability criterion
        print('QUICK stable')
    else:
        #Otherwise unstable
        print('unstable')

# name = 'FTCS'
# print(name)
# C = [0.05, 0.50, 0.40, 0.35, 0.5]
# s = [0.20, 0.25, 0.25, 0.40, 0.5]

# dx, dt = np.zeros(len(C)), np.zeros(len(C))
# for i, (Ci, si) in enumerate(zip(C, s)):
#     dx[i], dt[i] = Cs2dxdt(Ci,si)
# for dxi, dti, Ci, si in zip(dx, dt, C, s):
#         StabilityCriteria(name, Ci, si, dxi, dti)

# name = 'FTCS'
# print(name)
dx = [0.005, 0.02, 0.025, 0.04, 0.049]
dt = [0.002, 0.04, 0.0625, 0.08, 0.125]

C, s = np.zeros(len(dx)), np.zeros(len(dx))
for i, (dxi, dti) in enumerate(zip(dx, dt)):
    C[i], s[i] = dxdt2Cs(dxi,dti)

methods = ['FTCS', 'Upwind', 'Trapezoidal', 'QUICK']
for name in methods:
    print(name)
    for dxi, dti, Ci, si in zip(dx, dt, C, s):
        StabilityCriteria(name, Ci, si, dxi, dti)
