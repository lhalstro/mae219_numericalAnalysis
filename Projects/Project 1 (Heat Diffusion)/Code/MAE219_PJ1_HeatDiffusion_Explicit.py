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

#Plot Temperature Distribution
plt.figure()
plt.axis([0,L,T0,T1])
plt.title('Slab Temperature Distribution Over Time', fontsize=txt_ttl)
plt.xlabel('Length [nd]', fontsize=txt_lbl)
plt.ylabel('Temperature [nd]', fontsize=txt_lbl)

time = 0. #Initialize time counter
Told = Tnew #Initialize old temperature vector
while time <= t_end:

    for i in range(1,n-1):
        Tnew[i]= s*Told[i + 1] + (1-2.0*s)*Told[i] + s*Told[i - 1]

    plt.plot(x,Tnew,linewidth=1)
    time = time + dt
    Told = Tnew


plt.show()



##IMPLICIT SOLUTION#########

def trisolv(a,b,c,z):
    m=b.size
    for i in range(m-1):
        norm = a[i+1] / b[i]            #scale row to get zero result from subtraction
        a[i+1] = a[i+1] - norm * b[i]   #make below value zero
        b[i+1] = b[i+1] - norm * c[i]   #Complete row subtraction
        z[i+1] = z[i+1] - norm * z[i]   #Mirror subtraction in auxillary matrix

    #Solve last row
    x[m,1] = z[m] / b[m]
    #Zip up solution
    for i in range(m-1, 0, -1):
        x[i,1] = (z[i] - c[i] * x[i+1]) / b[i]

    return x


