#CASE STUDY 1 - 1D TRANSIENT HEAT TRANSFER
#Logan Halstrom
#MAE 219
#Created:  10/14/14

#Description:  Forward-Time, Centered-Space (FTCS) numeric solution of 1D unsteady heat transfer

import numpy as np
import matplotlib.pyplot as plt
import os

def SolvAna(x, t_end):
    #Heat diffusion analytic solution at all given points and a given time 
    #is an infinite sum.  Add terms to sum until residual is sufficiently small
    #Returns analytic solution for temperature distribution at a given time.

    n = len(x)
    Tana = [0]*n
    for i in range(0, n):
        #analytic solution
        ana = 1.
        for k in range(1, 100):
            kterm = 2. * k - 1.
            ana = ana - ((4. / (kterm * np.pi)) * 
                        np.sin(kterm * np.pi * x[i]) * 
                        np.exp( -(kterm) ** 2. * np.pi **2. * t_end))

            Tana[i] = ana


#        Tana[i] = ana

    return Tana


def SolvExp(Tnew, s, dt, t_end):
    #FTCS explicit scheme solution

    n = len(Tnew)
    Told = list(Tnew) #Initialize old temperature vector
    time = 0. #Initialize time counter
    while time <= t_end:

        for i in range(1, n-1):
            Tnew[i]= s*Told[i + 1] + (1-2.0*s)*Told[i] + s*Told[i - 1]

        time = time + dt
        Told = list(Tnew)
    return Tnew

def SolvImp(Tnew, s, dt, t_end):
    #FTCS implicit scheme solution

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

    Told = list(Tnew)
    time = 0. #Initialize time counter
    while time <= t_end:
        Tnew = trisolv(aa, bb, cc, Told)    #solve system at current time
        Told = list(Tnew)
        time += dt          #step forward in time
    return Tnew

def trisolv(a,b,c,d):
    #Tridiagonal Matrix Algorithm

   n = len(b)
   #initialize
   cprime, dprime, x = [0]*n, [0]*n, [0]*n
   cprime[0] = c[0] / b [0]
   dprime[0] = d[0] / b[0]
   for i in range(1, n-1):
       cprime[i] = c[i] / (b[i] - a[i] * cprime[i-1])
       dprime[i] = (d[i] - a[i] * dprime[i-1]) / (b[i] - a[i] * cprime[i-1])
   dprime[-1] = (d[-1] - a[-1] * dprime[-2]) / (b[-1] - a[-1] * cprime[-2])
   #Zip-up solution
   x[-1] = dprime[-1]
   for i in range(n-2, -1, -1):
       x[i] = dprime[i] - cprime[i] * x[i+1]
   return x

#    n=len(b)
#    for i in range(n-2):
#        norm = a[i+1] / b[i]            #scale row to get zero result from subtraction
#        a[i+1] = a[i+1] - norm * b[i]   #make below value zero
#        b[i+1] = b[i+1] - norm * c[i]   #Complete row subtraction
#        d[i+1] = d[i+1] - norm * d[i]   #Mirror subtraction in auxillary mat#rix
#    x = [0]*n
#    #Solve last row
#    x[-1] = d[-1] / b#[-1]
#    #Zip up solution
#    for i in range(n-2, -1#, -1):
#        x[i] = (d[i] - c[i] * x[i+1]) / b#[i]
#    return x

def RMSerror(num, ana):
    #Find RMS error of a numeric solution compared to the analytic solution

    n = len(num)
    rms = 0
    for i in range(0,n):
        rms += ( (num[i] - ana[i]) ** 2. )

    rms = (rms / n) ** (1. / 2.)

    return rms

def main():
    #Initialize variables and call solvers.  Plot RMS error comparison.

    #FIGURE SETUP
    WIDTH = 495.0  # the number latex spits out
    FACTOR = 1.0   # the fraction of the width you'd like the figure to occupy
    fig_width_pt  = WIDTH * FACTOR

    inches_per_pt = 1.0 / 72.27
    golden_ratio  = (np.sqrt(5) - 1.0) / 2.0  # because it looks good

    fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
    fig_height_in = fig_width_in * golden_ratio   # figure height in inches
    fig_dims      = [fig_width_in, fig_height_in] # fig dims as a list

    #INPUTS
    #Parametric study inputs:
    t_end = [0.03, 0.06, 0.09]  
    s = [1. / 6., 0.25, 0.5, 0.75]
    plot_temperature = 0        #plot temperature distribution figures = 1
    plot_rms = 0        #plot RMS error = 1
#    s = [1. / 6., .25]
#    t_end = [0.03, 0.06]
    #Other Inputs
    L = 1.              #Domain lenghth       [n.d.]
    T0 = 0.             #Initial temperature  [n.d.]
    T1 = 1.             #Boundary temperature [n.d.]
    t_start = 0.
    #s = 1. / 6          #dt/(dx)^2
    n = 21

    txt_lbl = 14                        #label fontsize
    txt_ttl = 14                        #title fontsize
    mkr = 8

    #MESH
    #dx        = L / (N - 1)
    x = np.linspace(0,L,n)
    dx = x[1]-x[0]
     
    #CALCULATE TIME STEP (Based on s and dx)
    #dt = np.power(np.multiply(s, dx), 2.)   

    #Initial Temperature Distribution
    Tinitial = [T0]*n
    #Boundary conditions
    Tinitial[0]  = T1
    Tinitial[-1] = T1

    #SOLUTIONS
    #Loop through input parameters
    size = (len(s), len(t_end))
    rms_exp, rms_imp = np.zeros(size), np.zeros(size)
    for i in range(0, len(s)):
#        rms_current = [0]*len(t_end)
        for j in range(0, len(t_end)):
            dt = s[i] * dx ** 2.
            #Analytic Solution
            Tanalytic = SolvAna(x, t_end[j])
            #Explicit Numerical Solution
            Texplicit = SolvExp(np.array(Tinitial), s[i], dt, t_end[j])
            #Implicit Numerical Solution
            Timplicit = SolvImp(np.array(Tinitial), s[i], dt, t_end[j])
            #Root Mean Square Error
            rms_exp[i,j] = RMSerror(Texplicit, Tanalytic)
            rms_imp[i,j] = RMSerror(Timplicit, Tanalytic)
            
            #PLOTTING
            plt.figure(figsize=fig_dims)
            #plt.axis([0,L,T0,T1])
            plt.title('Temperature Distribution at t=' + str(t_end[j]) + ' for s=' + str(s[i])[:5], fontsize=txt_ttl)
            plt.xlabel('Length [nd]', fontsize=txt_lbl)
            plt.ylabel('Temperature [nd]', fontsize=txt_lbl)
            plt.plot(x, Tanalytic, '-go', label='Analytic Solution', markersize=mkr+2, mfc='none', linewidth=1)
            plt.plot(x, Texplicit, '-rs', label='Explicit Solution', markersize=mkr, linewidth=1)
            plt.plot(x, Timplicit, '-b^', label='Implicit Solution', markersize=mkr, linewidth=1)
            plt.legend(loc='best')
            #Save
            save_name = 'PJ1_s' + str(s[i])[:5] + '_t' + str(t_end[j]) + '.pdf'
            try:
                os.mkdir('Results')
            except Exception:
                pass
            plt.savefig('Results/' + save_name, bbox_inches='tight') 
            
            if plot_temperature == 1:
                plt.show()
            plt.clf()

    #RMS ERROR
    #Save
    np.savetxt('Results/PJ1_RMSerror_Explicit.dat', rms_exp, delimiter=',')
    np.savetxt('Results/PJ1_RMSerror_Implicit.dat', rms_imp, delimiter=',')
    #PLOT RMS ERROR
    if plot_rms == 1:
        for i in range(0, len(t_end)):
            plt.figure()
            #plt.axis([0,L,T0,T1])
            plt.title('RMS Error at t = %s'%(t_end[i]), fontsize=txt_ttl)
            plt.xlabel('s [dt/dx^2]', fontsize=txt_lbl)
            plt.ylabel('RMS Error', fontsize=txt_lbl)
            plt.plot(s, rms_exp[:,i], 'r.', label='Explicit', linewidth=1)
            plt.plot(s, rms_imp[:,i], 'b.', label='Implicit', linewidth=1)
            plt.legend(loc='best')
            plt.show()
    

if __name__ == "__main__":
   main()