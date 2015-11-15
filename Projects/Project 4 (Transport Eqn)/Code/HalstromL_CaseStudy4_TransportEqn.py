"""
CASE STUDY 4 - TRANSPORT EQUATION
MAE 219
LOGAN HALSTROM
15 NOVEMBER 2014

DESCRIPTION:  This code uses a number of numerical schemes to solve the
transport equation.
"""

import numpy as np
import matplotlib.pyplot as plt 
import os
import time
import pandas as pd
from scipy import stats

# Configure figures for production
WIDTH = 495.0  # width of one column
FACTOR = 1.0   # the fraction of the width the figure should occupy
fig_width_pt  = WIDTH * FACTOR

inches_per_pt = 1.0 / 72.27
golden_ratio  = (np.sqrt(5) - 1.0) / 2.0      # because it looks good
fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
fig_height_in = fig_width_in * golden_ratio   # figure height in inches
fig_dims      = [fig_width_in, fig_height_in] # fig dims as a list

txt_tit = 18
txt_lbl = txt_tit
#Line Styles
mark = 8
line = 2.5
#Font Styles
font_tit = {'family' : 'serif',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 18,
            }
font_lbl = {'family' : 'serif',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 18,
            }
font_box = {'family' : 'arial',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 12,
            }
#Textbox Properties
textbox_props = dict(boxstyle='round', facecolor='white', alpha=0.5)
#Figure Save Filetype
filetype = 'pdf'

class SolnData:
    def __init__(self, name):
        # ni, nj = len(cases), len(x)
        self.name = name
        self.x = []
        self.phi = []
        self.stable = []  #stability criterion (stable if 1). Assessed later
        self.RMS = []       #RMS error of method for each case

def TexTable(filename, A, rows, cols, decimal_points=''):
    """Given matrix of data, column/row titles, write table to .tex file in
    LaTeX format.
    NOTES:  use formatters to put same text in each column entry (i.e. units)

    filename --> name of savefile
    A --> data to tablulate
    rows --> row titles
    cols --> column titles
    decimal_points --> number of decimal points in table entries (Default is given format)
    """
    nx, ny = A.shape

    #SEPARATE EACH COLUMN AND EACH ROW WITH A LINE
    lines=0
    col_sep = ' | ' if lines==1 else ' '
    #BOLD TITLES
    for i, r in enumerate(rows): rows[i] = '\\textbf{' + r + '}'
    for i, c in enumerate(cols): cols[i] = '\\textbf{' + c + '}'
    #BLANK LEFT COLUMN OPTION
    if len(cols)==ny:
        #If user did not provide a column title for the leftmost column:
        #insert blank column title for leftmost column
        cols = np.insert(cols, 0, '{}')
    
    with open(filename, 'w') as f:
        f.write('\\begin{tabular}{|c | ' + col_sep.join(['c'] * (len(cols)-1)) + '|}\n')
        f.write('\hline\n')
        f.write(' & '.join([str(col) for col in cols]) + ' \\\\\n')
        f.write('\hline\n')
        for i, row in enumerate(rows):
            X = []
            for x in A[i,:]:
                if x > 1e2 or x< 10**-(decimal_points-1):
                    #show value in scientific notation if it is too large or too small
                    fmt = '{:.' + str(decimal_points) + 'e}'
                else:
                    #show value in floating point with specified decimals
                    fmt = '{:.' + str(decimal_points) + 'f}'
                X.append(fmt.format(x))
            f.write(row + ' & ' + ' & '.join(X) + ' \\\\\n')
            if lines==1: f.write('\hline\n')        
        if lines !=1: f.write('\hline\n')
        f.write('\end{tabular}')

    # #Attempt at DataFrame version
    # df = pd.DataFrame(A, index=rows, columns=cols)
    
    # with open(filename, 'w') as f:
    #     f.write('\\begin{tabular}{' + ' | '.join(['c|'] * len(df.columns)-1) + '}\n')
    #     f.write('\hline\n')
    #     f.write(' & '.join([str(x) for x in df.columns.values]) + ' \\\\\n')
    #     for i, row in df.iterrows():
    #         if i==1:
    #             f.write('\hline\n')
    #         f.write(str(df.index.values[i]) + ' & ' + ' & '.join([str(x) for x in row.values]) + ' \\\\\n')
    #     f.write('\hline\n')
    #     f.write('\end{tabular}')

def MakeOutputDir(Dir):
    """make results output directory if it does not already exist.  Input directory name"""
    try:
        os.mkdir(Dir)
    except Exception:
        pass

def RMSerror(num, ana):
    """Find RMS error of a numeric solution compared to the analytic solution"""
    n = len(num)
    rms = 0
    for i in range(0,n):
        rms += ( (num[i] - ana[i]) ** 2. )
    rms = (rms / n) ** (1. / 2.)
    return rms

def TriPerSolv(a, b, c, d):
    """Solve periodic tridiagonal systems with cyclic Thomas Algorithm
    and Sherman-Morrison formula.
    """
    def TriSolv(a, b, c, d):
    #Tridiagonal Matrix Algorithm

       n = len(b)
       #initialize
       cprime, dprime, x = [0]*n, [0]*n, [0]*n
       #First Row
       cprime[0] = c[0] / b [0]
       dprime[0] = d[0] / b[0]
       #Middle Rows
       for i in range(1, n-1):
           cprime[i] = c[i] / (b[i] - a[i] * cprime[i-1])
           dprime[i] = (d[i] - a[i] * dprime[i-1]) / (b[i] - a[i] * cprime[i-1])
       #Last Row
       dprime[-1] = (d[-1] - a[-1] * dprime[-2]) / (b[-1] - a[-1] * cprime[-2])
       #Zip-up solution
       x[-1] = dprime[-1]
       for i in range(n-2, -1, -1):
           x[i] = dprime[i] - cprime[i] * x[i+1]
       return x

    n = len(d)
    #SPLIT TRIDIAGONAL MATRIX FROM PERIODIC DRIDIAGONAL MATRIX
    aa = [a] * n
    bb = [b] * n
    cc = [c] * n
    #Boundary Conditions
    aa[0], bb[0], cc[0]= 0, b-a, c
    aa[-1], bb[-1], cc[-1]= a, b-c, 0
    #OTHER HALF OF SPLIT
    uu = np.zeros(n)
    uu[0], uu[-1] = a, c
    vv = np.zeros(n)
    vv[0], vv[-1] = 1, 1

    #SOLVE TRIDIAGONAL SYSTEMS
    z = TriSolv(aa, bb, cc, uu)
    y = TriSolv(aa, bb, cc, d)

    x = np.subtract(y, 1/(1 + np.dot(vv, z)) 
        * np.dot(z, np.dot(vv, y)))
    return x


def AnaSoln(x, t):
    """analytic solution
    x --> 1D space mesh vector
    t --> single time to evaluate at
    """
    phi = np.exp(-k**2 * D * t) * np.sin(k * np.subtract(x, u*t))
    return phi

def InitialCondition(x):
    """ Initial condition (Phi(x,0)) for numeric solutions"""
    return np.sin(k * x)

def FTCS(Fini, tend, dt, s, C):
    """Forward Time Central Space explicit solution of the transport eqn.
    Fini --> initial condition for Phi
    tend --> end time for simulation
    """
    n = len(Fini)
    Fold = np.array(Fini)
    Fnew = np.empty_like(Fold)
    t = 0
    while t < tend:
        #Periodic BCs
        # Fnew[0] = (s * (Fold[1] - 2*Fold[0] + Fold[-1]) + Fold[0]
        #                 -C/2 * (Fold[1] - Fold[-1]))
        for i in range(1, n-1):
            Fnew[i] = (s * (Fold[i+1] - 2*Fold[i] + Fold[i-1]) + Fold[i]
                        -C/2 * (Fold[i+1] - Fold[i-1]))
        Fnew[-1] = (s * (Fold[1] - 2*Fold[-1] + Fold[-2]) + Fold[-1]
                        -C/2 * (Fold[1] - Fold[-2]))
        Fnew[0] = Fnew[-1]
        t += dt
        Fold = np.array(Fnew)
    return Fnew

def Upwind(Fini, tend, dt, s, C):
    """Explicit forward Euler solution of the transport eqn, 
    with second-order upwind scheme for convective flux and 
    central differencing for diffusive flux.
    Fini --> initial condition for Phi
    tend --> end time for simulation
    """
    n = len(Fini)
    Fold = np.array(Fini)
    Fnew = np.empty_like(Fold)
    t = 0
    while t < tend:
        #backward difference for u>0
        if u>0:
            #Periodic BCs
            # Fnew[0] = (s * (Fold[1] - 2*Fold[0] + Fold[-1]) + Fold[0]
            #                 -C/2 * (3*Fold[0]  - 4*Fold[-1] + Fold[-2]))
            Fnew[1] = (s * (Fold[2] - 2*Fold[1] + Fold[-1]) + Fold[1]
                            -C/2 * (3*Fold[1]  - 4*Fold[-1] + Fold[-2]))
            for i in range(2, n-1):
                Fnew[i] = (s * (Fold[i+1] - 2*Fold[i] + Fold[i-1]) + Fold[i]
                            -C/2 * (3*Fold[i]  - 4*Fold[i-1] + Fold[i-2]))
            Fnew[-1] = (s * (Fold[1] - 2*Fold[-1] + Fold[-2]) + Fold[-1]
                            -C/2 * (3*Fold[-1]  - 4*Fold[-2] + Fold[-3]))
            Fnew[0] = Fnew[-1]
        
        #forward difference for u<0
        elif u<0:
            #Periodic BCs
            # Fnew[0] = (s * (Fold[1] - 2*Fold[0] + Fold[-1]) + Fold[0]
            #                 -C/2 * (-Fold[2] + 4*Fold[1] - 3*Fold[0]))
            for i in range(1, n-2):
                Fnew[i] = (s * (Fold[i+1] - 2*Fold[i] + Fold[i-1]) + Fold[i]
                            -C/2 * (-Fold[i+2] + 4*Fold[i+1] - 3*Fold[i]))
            Fnew[-2] = (s * (Fold[-1] - 2*Fold[-2] + Fold[-3]) + Fold[-2]
                            -C/2 * (-Fold[1] + 4*Fold[-1] - 3*Fold[-2]))
            Fnew[-1] = (s * (Fold[0] - 2*Fold[-1] + Fold[-2]) + Fold[-1]
                            -C/2 * (-Fold[2] + 4*Fold[1] - 3*Fold[-1]))
            Fnew[0] = Fnew[-1]

        t += dt
        Fold = np.array(Fnew)
    return Fnew

def CrankNicholson(Fini, tend, dt, s, C, theta=0.5):
    """Trapezoidal Method (Crank-Nicholson). Weighted average of 
    explicit and implicit methods.
    Fini --> initial condition for Phi
    tend --> end time for simulation
    theta --> weight of each method (default 1/2 for true average)
    """
    #Use Cyclic Thomas Tridiagonal Matrix Algorithm (1) or linalg.solve (0)
    CTDMA = 1

    n = len(Fini)
    #Coefficients for periodic tridiagonal matrix
    aa = theta*(C-s)
    bb = 1 + theta*2*s
    cc = theta*(-s-C)

    current_method = 'CTDMA'
    if CTDMA != 1:
        current_method = 'Gaussian Elimination'
        #Build periodic matrix for linalg.solve
        size = (n,n)
        A = np.zeros(size)
        A[0,0], A[0,1], A[0,-1] = bb, cc, aa
        A[-1,0], A[-1,-2], A[-1,-1] = cc, aa, bb
        for i in range(1, n-1):
            A[i,i-1], A[i,i], A[i,i+1] = aa, bb, cc
        # print(A)
    
    Fold = np.array(Fini)
    Fnew = np.empty_like(Fold)
    t = 0
    #measure wall clock time of solution to compare Gaussian Elimination with CTDMA
    timestart = time.time()
    while t < tend:
        #RHS
        d = np.zeros(n)
        d[0] = ((1-theta) * ((-C+s)*Fold[1] 
                    + (1/(1-theta)-2*s)*Fold[0] + (s+C)*Fold[-1]))
        for i in range(1, n-1):
            d[i] = ((1-theta) * ((-C+s)*Fold[i+1] 
                    + (1/(1-theta)-2*s)*Fold[i] + (s+C)*Fold[i-1]))
        d[-1] = ((1-theta) * ((-C+s)*Fold[0] 
                    + (1/(1-theta)-2*s)*Fold[-1] + (s+C)*Fold[-2]))

        if CTDMA == 1:
            #Use Cyclic Thomas Tridiagonal Matrix Algorithm
            Fnew = TriPerSolv(aa, bb, cc, np.array(d))
        else:
            #Use gaussian elimination
            Fnew = np.linalg.solve(A,d)

        t += dt
        Fold = np.array(Fnew)
    timeend = time.time()
    timerun = timeend - timestart
    print(current_method + ' Wall Clock Time: %1.4f'%(timerun))
    return Fnew

def QUICK(Fini, tend, dt, s, C):
    """Quadratic Upwind Interpolation for Convection Kinematics (QUICK) scheme.
    2nd order, non-oscillatory method.  Explicit forward euler in time,
    QUICK for convective flux (positive direction only), 
    central difference for diffusive flux.
    Fini --> initial condition for Phi
    tend --> end time for simulation
    """
    n = len(Fini)
    Fold = np.array(Fini)
    Fnew = np.empty_like(Fold)
    t = 0
    while t < tend:
        #Positive direction QUICK scheme only
        # Fnew[0] = (s*(Fold[1] - 2*Fold[0] + Fold[-1]) + Fold[0]
        #                 -C*(3/8*Fold[1] + 3/8*Fold[0] - 7/8*Fold[-1] + 1/8*Fold[-2]))
        Fnew[1] = (s*(Fold[2] - 2*Fold[1] + Fold[-1]) + Fold[1]
                        -C*(3/8*Fold[2] + 3/8*Fold[1] - 7/8*Fold[-1] + 1/8*Fold[-1]))
        for i in range(2, n-1):
            Fnew[i] = (s*(Fold[i+1] - 2*Fold[i] + Fold[i-1]) + Fold[i]
                        -C*(3/8*Fold[i+1] + 3/8*Fold[i] - 7/8*Fold[i-1] + 1/8*Fold[i-2]))
        Fnew[-1] = (s*(Fold[1] - 2*Fold[-1] + Fold[-2]) + Fold[-1]
                        -C*(3/8*Fold[1] + 3/8*Fold[-1] - 7/8*Fold[-2] + 1/8*Fold[-3]))
        Fnew[0] = Fnew[-1]

        t += dt
        Fold = np.array(Fnew)
    return Fnew

def Solvers(f, x, Fini, tend, dt, s, C):
    """input solution method object and case conditions,
    solve case using given method, and save result to object.
    """
    if f.name=='FTCS':
        f.x.append(x)
        f.phi.append(FTCS(np.array(Fini), tend, dt, s, C))
    elif f.name=='Upwind':
        f.x.append(x)
        f.phi.append(Upwind(np.array(Fini), tend, dt, s, C))
    elif f.name=='Trapezoidal':
        f.x.append(x)
        f.phi.append(CrankNicholson(np.array(Fini), tend, dt, s, C))
    elif f.name=='QUICK':
        f.x.append(x)
        f.phi.append(QUICK(np.array(Fini), tend, dt, s, C))
    return f

def StabilityCriteria(f, C, s, dx, dt):
    """Assess stability of given numerical method for given case
    f --> array of data for given method
    dx --> grid spacing for given case
    dt --> time spacing for given case
    """
    #SELECT STABILITY CRITERIA FOR SPECIFIC METHOD
    if f.name=='FTCS' and dx<=2*D/u and dt<=dx**2/(2*D):
        #FTCS stability criterion
        f.stable.append(1)
    elif f.name=='Upwind' and C+2*s<1:
        #Upwind stability criterion (CFL)
            f.stable.append(1)
    elif f.name=='Trapezoidal':
        #Trapezoidal is implicit, thus inhernetly stable
        f.stable.append(1)
    elif f.name=='QUICK' and C<=min(2-4*s, (2*s)**0.5):
        #QUICK stability criterion
        f.stable.append(1)
    else:
        #Otherwise unstable
        f.stable.append(0)
    return f

def PlotOneMethod(F, cases, t, showplot):
    """plot all cases for each method.
    F --> array of data for each method
    """
    n = len(F)-1  #number of plots is one less than length because of analytic soln
    #PLOT EACH METHOD
    print('')
    for imeth in range(0, n):
        print('Plotting Method %s of %s'%(imeth+1, n))

        #PLOT ALL CASES UNBOUND AND ONLY STABLE CASES BOUNDED
        for k in range(0, 2):
            #first loop bound by numeric solution, second loop bound by analytic
            fig, ax = plt.subplots(1)
            ax.set_title('1-D Transport Equation ' + str(F[imeth].name) + ' Solution (t=%1.2fs)'%(t), fontdict=font_tit)
            ax.set_xlabel(r'$x\,[m]$', fontdict=font_lbl)
            ax.set_ylabel(r'$\Phi$', fontdict=font_lbl)
            ax.plot(F[-1].x[0], F[-1].phi[0], '--', label='Analytic', linewidth=line, markersize=mark)
            # ymin = ymax = 0
            unstable = 'Unstable Cases:'
            for i, (C,s) in enumerate(cases):
                # #Find y-bounds of numeric solution within given range
                # c1 = 1e2
                # lobound, upbound = c1*F[-1].phi[0].min(), c1*F[-1].phi[0].max()
                # curmin, curmax = F[imeth].phi[i].min(), F[imeth].phi[i].max()
                # if curmin<ymin and curmin>lobound: ymin = curmin
                # if curmax>ymax and curmax<upbound: ymax = curmax
                # if (curmin<lobound or curmax>upbound) and k==1:
                #     #if solution exceeds bounds and this is bounded plot, don't plot
                #     continue
                if F[imeth].stable[i]!=1:
                    #Add case to list of unstable cases
                    unstable += '\nC=%1.2f, s=%1.2f'%(C, s)
                    if k==1:
                        #don't plot case if it is unstable and if this is the bounded plot
                        continue
                ax.plot(F[imeth].x[i], F[imeth].phi[i], '-.', 
                    label='C=%1.2f, s=%1.2f'%(C, s), 
                    linewidth=line, markersize=mark)
            # if k==1:
            #     #bound axis proportinally to analytic solution
            #     c2 = 1
            #     plt.ylim(c2*ymin, c2*ymax)
            ax.legend(loc='upper right', fancybox=True, framealpha=0.5)
            textboxtext = unstable
            ax.text(0.05, 0.05, textboxtext, transform=ax.transAxes, 
                verticalalignment='bottom', horizontalalignment='left',
                fontdict=font_box, bbox=textbox_props)
            MakeOutputDir('Results')
            if k==0:
                save_name = str(F[imeth].name) + '_t%1.2f.'%(t) + filetype
            elif k==1:
                save_name = str(F[imeth].name) + '_t%1.2f_AxLim.'%(t) + filetype
            plt.savefig('Results/' + save_name, bbox_inches='tight')
            if showplot!=0:
                plt.show()
            plt.close()

def PlotAllMethods(F, cases, t, showplot):
    """plot all cases for one method.
    F --> array of data for each method
    imeth --> index of method to plot
    """
    print('')
    for i, (C, s) in enumerate(cases):
        #PLOT TRANSPORT EQUATION
        print('Plotting Case %s of %s'%(i+1, len(cases)))
        plt.figure()
        plt.title('1-D Transport Equation (t=%1.2fs, C=%1.2f, s=%1.2f)'%(t, C, s), fontdict=font_tit)
        plt.xlabel(r'$x\,[m]$', fontdict=font_lbl)
        plt.ylabel(r'$\Phi$', fontdict=font_lbl)
        plt.plot(F[-1].x[0], F[-1].phi[0], '--', label='Analytic', linewidth=line, markersize=mark)
        for j in range(0, len(F)-1):
            plt.plot(F[j].x[i], F[j].phi[i], '.', label=F[j].name, 
                linewidth=line, markersize=mark)

        plt.legend(loc='best', fancybox=True, framealpha=0.5)
        MakeOutputDir('Results')
        save_name = 'AllCases_t%1.2f_C%1.2f_s%1.2f.'%(t, C, s) + filetype
        plt.savefig('Results/' + save_name, bbox_inches='tight')
        if showplot==1:
            plt.show()
        plt.close()

        #PLOT ERROR COMPARED TO ANALYTIC SOLUTION
        # print('Plotting Error')
        plt.figure()
        plt.title('Numerical Error (t=%1.2fs, C=%1.2f, s=%1.2f)'%(t, C, s), fontdict=font_tit)
        plt.xlabel(r'$x\,[m]$', fontdict=font_lbl)
        plt.ylabel(r'$\Phi_{ana}-\Phi_{num}$', fontdict=font_lbl)
        for j in range(0, len(F)-1):
            plt.plot(F[j].x[i], F[j].phi[i]-F[-1].phi[i], '.', label=F[j].name, 
                linewidth=line, markersize=mark)

        plt.legend(loc='best', fancybox=True, framealpha=0.5)
        MakeOutputDir('Results')
        save_name = 'AllCases_Error_t%1.2f_C%1.2f_s%1.2f.'%(t, C, s) + filetype
        plt.savefig('Results/' + save_name, bbox_inches='tight')
        if showplot==1:
            plt.show()
        plt.close()

def PlotGridSensitivity(F, cases, DX, DT, t, showplot):
    """for sets of constant dx or dt, plot error to demonstrate sensitivity to grid size
    or time step.  Plot log error to show order accuracy.
    F --> array of data for each method
    """
    #ORGANIZE DATA SETS
    #Most common value of dx/dt (mode)
    dxconst, dtconst = stats.mode(DX)[0][0], stats.mode(DT)[0][0]
    idt, idx = np.where(DX==dxconst)[0], np.where(DT==dxconst)[0]   #indicies to plot
    print('')
    #CONSTANT DX CASE, PLOT AGAINST DT
    if len(idt)>2:
        print('Plotting Time Sensitivity')
        # idt = np.where(DX==dxconst)[0]  #indicies to plot
        plt.figure()
        plt.title('Time Step Sensitivity (t=%1.2fs)'%(t), fontdict=font_tit)
        plt.xlabel(r'$\Delta t\,[s]$', fontdict=font_lbl)
        plt.ylabel(r'$RMS_{error}$', fontdict=font_lbl)
        for j in range(0, len(F)-1):
            plt.plot([DT[i] for i in idt], [F[j].RMS[i] for i in idt], 
                '.', label=F[j].name, linewidth=line, markersize=mark)

        plt.legend(loc='best', fancybox=True, framealpha=0.5)
        plt.ylim([0,10])
        MakeOutputDir('Results')
        save_name = 'TimeSens_t%1.2f_dt%1.4f-%1.2f_n%1.0f.'%(t, DT[idt[0]], DT[idt[-1]], len(idt)) + filetype
        plt.savefig('Results/' + save_name, bbox_inches='tight')
        if showplot==1:
            plt.show()

        for f in F:
            if f.name=='Analytic': continue
            plt.figure()
            xplot, yplot = [DT[i] for i in idt], [f.RMS[i] for i in idt]
            #slope of log-log plot is order of accuracy
            slope, intercept = np.polyfit(np.log(xplot), np.log(yplot), 1)
            plt.title(f.name + ' Effective Order of Accuracy in Time: %1.2f'%(slope), fontdict=font_tit)
            plt.xlabel(r'$\Delta t\,[s]$', fontdict=font_lbl)
            plt.ylabel(r'$RMS_{error}$', fontdict=font_lbl)
            plt.loglog(xplot, yplot, '.', label=f.name, linewidth=line, markersize=mark)

            # plt.legend(loc='best', fancybox=True, framealpha=0.5)
            # plt.ylim([0,10])
            MakeOutputDir('Results')
            save_name = 'TimeSens_Log_' + f.name + '_t%1.2f_dt%1.4f-%1.2f_n%1.0f.'%(t, DT[idt[0]], DT[idt[-1]], len(idt)) + filetype
            plt.savefig('Results/' + save_name, bbox_inches='tight')
            if showplot==1:
                plt.show()

    #CONSTANT DT CASE, PLOT AGAINST DX
    if len(idx)>2:
        print('Plotting Mesh Sensitivity')
        #Constant dx case, plot against dt
        # idx = np.where(DT==dtconst)[0]  #indicies to plot
        plt.figure()
        plt.title('Mesh Sensitivity (t=%1.2fs)'%(t), fontdict=font_tit)
        plt.xlabel(r'$\Delta x\,[m]$', fontdict=font_lbl)
        plt.ylabel(r'$RMS_{error}$', fontdict=font_lbl)
        for j in range(0, len(F)-1):
            plt.plot([DX[i] for i in idx], [F[j].RMS[i] for i in idx], 
                '.', label=F[j].name, linewidth=line, markersize=mark)

        plt.legend(loc='best', fancybox=True, framealpha=0.5)
        MakeOutputDir('Results')
        save_name = 'MeshSens_t%1.2f_dt%1.4f-%1.2f_n%1.0f.'%(t, DX[idt[0]], DX[idt[-1]], len(idx)) + filetype
        plt.savefig('Results/' + save_name, bbox_inches='tight')
        if showplot==1:
            plt.show()

def PlotLog(F, cases, DX, DT, t, showplot):
    for f in F:
        if f.name=='Analytic': continue
        plt.figure()
        xplot, yplot = DT, f.RMS
        #slope of log-log plot is order of accuracy
        slope, intercept = np.polyfit(np.log(xplot), np.log(yplot), 1)
        plt.title(f.name + ' Effective Order of Accuracy in Time: %1.2f'%(slope), fontdict=font_tit)
        plt.xlabel(r'$\Delta t\,[s]$', fontdict=font_lbl)
        plt.ylabel(r'$RMS_{error}$', fontdict=font_lbl)
        plt.loglog(xplot, yplot, '.', label=f.name, linewidth=line, markersize=mark)

        # plt.legend(loc='best', fancybox=True, framealpha=0.5)
        # plt.ylim([0,10])
        MakeOutputDir('Results')
        save_name = 'TimeAcc_Log_' + f.name + '_t%1.2f_dt%1.4f-%1.2f_n%1.0f.'%(t, DT[0], DT[-1], len(DT)) + '.pdf'
        plt.savefig('Results/' + save_name, bbox_inches='tight')
        if showplot==1:
            plt.show()

        for f in F:
            if f.name=='Analytic': continue
            plt.figure()
            xplot, yplot = DX, f.RMS
            #slope of log-log plot is order of accuracy
            slope, intercept = np.polyfit(np.log(xplot), np.log(yplot), 1)
            plt.title(f.name + ' Effective Order of Accuracy in Space: %1.2f'%(slope), fontdict=font_tit)
            plt.xlabel(r'$\Delta x\,[s]$', fontdict=font_lbl)
            plt.ylabel(r'$RMS_{error}$', fontdict=font_lbl)
            plt.loglog(xplot, yplot, '.', label=f.name, linewidth=line, markersize=mark)

            # plt.legend(loc='best', fancybox=True, framealpha=0.5)
            # plt.ylim([0,10])
            MakeOutputDir('Results')
            save_name = 'MeshAcc_Log_' + f.name + '_t%1.2f_dt%1.4f-%1.2f_n%1.0f.'%(t, DX[0], DX[-1], len(DX)) + '.pdf'
            plt.savefig('Results/' + save_name, bbox_inches='tight')
            if showplot==1:
                plt.show()
            plt.close()

def PlotRMSActual(F, cases, t, showplot):
    
    #SAVE DATA IN LATEX TABLE FORMAT
    rows = []
    scols = []
    Ccols = []
    for C,s in cases:
        if s == 0.25:
            #Constant diffusion coefficient
            scols.append('%1.2f'%(C))
        if C == 0.5:
            Ccols.append('%1.2f'%(s))

    size = (len(F)-1,3)
    sTable, CTable = np.zeros(size), np.zeros(size)
    for i, f in enumerate(F):
        if f.name=='Analytic': continue #Leave out analytic solution
        #Method name for row title
        rows.append(f.name)
        k, l = 0, 0
        for j, (C,s) in enumerate(cases):
            if s == 0.25:
                #Constant diffusion coefficient
                # scols[k] = '%1.2f'%(C)
                sTable[i,k] = f.RMS[j]
                k += 1
                # print(scols)
            if C == 0.5:
                #constant courant number
                # Ccols.append('%1.2f'%(s))
                CTable[i,l] = f.RMS[j]
                l += 1
                # print(Ccols)
    #save tables
    TexTable('Results/RMS_sConst.tex', sTable, list(rows), scols, 3)
    TexTable('Results/RMS_CConst.tex', CTable, list(rows), Ccols, 3)

    #STABILITY TABLE
    cols = []
    size = (len(cases),len(F)-1)
    stabTable = np.zeros(size)
    RMStable = np.zeros(size)
    for j, (C,s) in enumerate(cases):
        cols.append('C=%1.2f, s=%1.2f'%(C,s))
        for i, f in enumerate(F):
            if f.name=='Analytic': continue
            stabTable[j,i] = f.stable[j]
            RMStable[j,i] = f.RMS[j]

    TexTable('Results/Stability.tex', stabTable, list(cols), list(rows), 4)
    TexTable('Results/RMS.tex', RMStable, list(cols), list(rows), 3)


    return

def main(cases, methods):
    """Main function, loop for various sets of inputs
    cases --> array of tuples of C=u*dt/dx and s=D*dt/dx^2
    method --> name of method for saving purposes
    actual --> 1 for case set given in problem statement.  Controls specific funcitons
    """
    global u, D, L, k, tend
    u = 0.2                 #convection velocity [m/s]
    D = 0.005               #diffusion coefficient [m^2/s]
    L = 1.                  #length [m]
    k = 2*np.pi/L           #constant [m^-1]
    tend = 1./(k**2 * D)     #tau (max time) [s]
    #TIME TO SOLVE FOR:
    t=tend

    #CHECK IF THIS IS SET OF CASES PRESCRIBED BY PROBLEM STATEMENT
    if (methods==['FTCS', 'Upwind', 'Trapezoidal', 'QUICK'] and
    cases==[(0.1, 0.25), (0.5, 0.25), (2., 0.25), (0.5, 0.5), (0.5, 1.)]):
        actual = 1
    else:
        actual = 0

    #INITIALIZE SOLUTIONS
    F = np.zeros(len(methods)+1, dtype=object)
    F[-1] = SolnData('Analytic')
    for i, method in enumerate(methods):
        F[i] = SolnData(method)
    DX, DT = np.zeros(len(cases)), np.zeros(len(cases))

    #SOLVE EACH CASE
    print('')
    for i, (C,s) in enumerate(cases):
        print('Solving case %s of %s'%(i+1, len(cases)))
        #MESH PARAMETERS
        dx = (D*C) / (u*s)
        x = np.append(np.arange(0, L, dx), L)
        nx = len(x)
        dt = s * dx**2 / D
        #store mesh spacing for each case
        DX[i], DT[i] = dx, dt

        #INITIAL CONDITION
        Fini = InitialCondition(x)

        #ANALYTIC SOLUTION
        F[-1].x.append(x)
        F[-1].phi.append(AnaSoln(x,t))
        #ALL NUMERIC SOLUTIONS
        for f in F:
            if f.name=='Analytic': continue #Leave out analytic solution
            #SOLVE
            f = Solvers(f, x, np.array(Fini), t, dt, s, C)
            #ASSESS STABILITY CRITERIA
            f = StabilityCriteria(f, C, s, dx, dt)
            #RMS ERROR
            f.RMS.append(RMSerror(f.phi[i], F[-1].phi[i]))

    #PLOT ALL CASES FOR EACH METHOD
    PlotOneMethod(F, cases, t, 0)
    #PLOT ALL METHODS FOR EACH CASE
    PlotAllMethods(F, cases, t, 0)

    if actual==1:
        #RMS ERROR FOR PRESCRIBED CASES
        print('\nAssessing RMS of prescribed case')
        PlotRMSActual(F, cases, t, 0)

    #PLOT MESH/TIME STENSITIVITY STUDY
    PlotGridSensitivity(F, cases, DX, DT, t, 0)
    PlotLog(F, cases, DX, DT, t, 0)

    print(DX, DT)

    print('\nDone!')

if __name__ == "__main__":

    #Actual Run
    methods = ['FTCS', 'Upwind', 'Trapezoidal', 'QUICK']
    cases = [(0.1, 0.25), (0.5, 0.25), (2., 0.25), (0.5, 0.5), (0.5, 1.)]

    # #Mesh sensitivity study
    # dx = [0.005, 0.03125, 0.125, 0.625, 2]
    # dt = [0.005] * len(dx)
    # #Time step sensitivity study
    # dt = [0.005, 0.03125, 0.125, 0.625, 2]
    # dt = [0.005, 0.03125, 0.125, 0.625, 2]
    # dx = [0.01] * len(dt)
    # Cs, ss = 0.2 * np.divide(dt,dx), 0.005 * np.divide(dt,np.power(dx,2))
    # cases = list(zip(Cs,ss))

    # C = [0.10, 0.50, 0.40, 0.35, 0.5]
    # s = [0.25, 0.25, 0.25, 0.40, 0.5]

    # C = [0.1, 0.2, 0.3, 0.05, 0.1]
    # s = [0.4, 0.3, 0.2,  0.1, 0.1]

    # C = [0.5, 0.6, 0.7, 0.8, 0.9]
    # s = [0.25, 0.25, 0.25, 0.25, 0.25]

    # C = [0.25, 0.4, 0.5, 0.6,  0.7]
    # s = [0.25, 0.25, 0.25, 0.25, 0.25]
    # cases = list(zip(C,s))


    main(cases, methods)