#CASE STUDY 1 - 1D TRANSIENT HEAT TRANSFER
#Logan Halstrom
#MAE 219
#Created:  10/14/14

#Description:  Forward-Time, Centered-Space (FTCS) numeric solution of 1D unsteady heat transfer

import numpy as np
import matplotlib.pyplot as plt

#INPUTS
L = 1.              #Domain lenghth       [n.d.]
T0 = 0.             #Initial temperature  [n.d.]
T1 = 1.             #Boundary temperature [n.d.]
t_start = 0.
t_end = 0.1
s = 1. / 6          #dt/(dx)^2
n = 21

txt_lbl = 14                        #label fontsize
txt_ttl = 14                        #title fontsize

#MESH
#dx        = L / (N - 1)
x = np.linspace(0,L,n)
dx = x[1]-x[0]
 
#CALCULATE TIME STEP (Based on s and dx)
dt = s*dx**2.0   

##EXPLICIT SOLVER####################

# Initial Condition
Tnew = [T0]*n

#Boundary conditions
Tnew[0]   = T1
Tnew[n-1] = T1


def ExplicitSoln(Tnew, t_end, dt, s):

    n = len(Tnew)
    Told = Tnew #Initialize old temperature vector
    time = 0. #Initialize time counter
    while time <= t_end:

        for i in range(1,n-1):
            Tnew[i]= s*Told[i + 1] + (1-2.0*s)*Told[i] + s*Told[i - 1]

        time = time + dt
        Told = Tnew
    return Tnew

def ImplicitSoln(Tnew, t_end, dt, s):
    
    n = len(Tnew)
    T1 = Tnew[0]        #Boundary temperature
    #Initialize tridiagonal matrix vectors
    aa = [-s] * n
    bb = [1 + 2*s] * n
    cc = [-s] * n
    #Boundary Conditions (b diagonal temperature known, a/c=0)
    aa[0], aa[-1] = 0, 0
    bb[0], bb[-1] = T1, T1
    cc[0], cc[-1] = 0, 0

    Told = Tnew
    time = 0. #Initialize time counter
    while time <= t_end:
        Tnew = trisolv(aa, bb, cc, Told)    #solve system at current time
        Told = Tnew
        time += dt          #step forward in time
    return Tnew


#Plot Temperature Distribution
plt.figure()
plt.axis([0,L,T0,T1])
plt.title('Slab Temperature Distribution Over Time', fontsize=txt_ttl)
plt.xlabel('Length [nd]', fontsize=txt_lbl)
plt.ylabel('Temperature [nd]', fontsize=txt_lbl)
plt.plot(x,Tnew,linewidth=1)
plt.show()



##IMPLICIT SOLUTION#########

def trisolv(a,b,c,d):
    n=b.size
    for i in range(n-1):
        norm = a[i+1] / b[i]            #scale row to get zero result from subtraction
        a[i+1] = a[i+1] - norm * b[i]   #make below value zero
        b[i+1] = b[i+1] - norm * c[i]   #Complete row subtraction
        d[i+1] = d[i+1] - norm * d[i]   #Mirror subtraction in auxillary matrix

    #Solve last row
    x[n,1] = d[n] / b[n]
    #Zip up solution
    for i in range(n-1, 0, -1):
        x[i,1] = (d[i] - c[i] * x[i+1]) / b[i]

    return x


