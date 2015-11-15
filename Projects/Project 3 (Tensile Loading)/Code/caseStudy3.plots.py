"""
CASE STUDY 3 - PERFORATED PLATE IN TENSION PLOTS
MAE 219
LOGAN HALSTROM
15 NOVEMBER 2014

DESCRIPTION:  This code plots results from a multiple runs for a plate 
under tension.
"""

import numpy as np
import matplotlib.pyplot as plt 
import os

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
mark = 8
line = 2

def SavePlot(name):
    try:
        os.mkdir('Results')
    except Exception:
        pass
    plt.savefig('Results/' + name)

def MakeGrid(R, L):
    n = 101
    r = np.linspace(R, L, n)
    thetamax = np.pi/2.
    theta = np.linspace(0, thetamax, n)
    rr, tt = np.meshgrid(r, theta)
    return r, theta, rr, tt

def VertSigmaXX(y):
    """Analytic solution for sigmaxx at x=0 symmetry plane
    sigma --> stress applied in x direction
    R --> radius of plate perforation
    y --> y axis vector
    """
    return sigma * (1 - 0.5*(R**2 / y**2) + 1.5*(R**4/y**4))

def VertSigmaYY(y):
    """Analytic solution for sigmaxx at x=0 symmetry plane
    sigma --> stress applied in x direction
    R --> radius of plate perforation
    y --> y axis vector
    """
    return sigma * (1.5*(R**2 / y**2) - 1.5*(R**4/y**4))

def VertSigmaXY(y):
    """Analytic solution for sigmaxx at x=0 symmetry plane
    sigma --> stress applied in x direction
    R --> radius of plate perforation
    y --> y axis vector
    """
    return sigma * (1.5*(R**2 / y**2) - 1.5*(R**4/y**4))

def HorzSigmaXX(x):
    """Analytic solution for sigmaxx at x=0 symmetry plane
    sigma --> stress applied in x direction
    R --> radius of plate perforation
    y --> y axis vector
    """
    return sigma * (1 - 2.5*(R**2 / x**2) + 1.5*(R**4/x**4))

def HorzSigmaYY(x):
    """Analytic solution for sigmaxx at x=0 symmetry plane
    sigma --> stress applied in x direction
    R --> radius of plate perforation
    y --> y axis vector
    """
    return sigma * (0.5*(R**2 / x**2) - 1.5*(R**4/x**4))

def DiagSigmaXX(x):
    """Analytic solution for sigmaxx at x=0 symmetry plane
    sigma --> stress applied in x direction
    R --> radius of plate perforation
    y --> y axis vector
    """
    return sigma * (1 + 0.5*(R**2 / x**2) - 3/8*(R**4/x**4))

def DiagSigmaYY(x):
    """Analytic solution for sigmaxx at x=0 symmetry plane
    sigma --> stress applied in x direction
    R --> radius of plate perforation
    y --> y axis vector
    """
    return sigma * (-0.5*(R**2 / x**2) + 3/8*(R**4/x**4))

class CaseData:
    def __init__(self, L, dx2):
        """Gathers data from cases and conditions
        for plotting."""
        self.L = L
        self.dx2 = dx2
        #Case Name
        self.case = str(L) + '_' + str(self.dx2) + '_case'
        self.label = 'FV Solution (L=' + str(self.L*2) + 'm, n$_{2}$=' + str(self.dx2) + ')'
        self.datapath = self.case + '/postProcessing/sets/100/'
        #Get Vertical Symmetry Plane Data
        raw = np.loadtxt(self.datapath + 'leftPatch_sigmaxx_sigmaxy_sigmayy.xy')
        self.y = raw[:,0]
        self.sXXvert = raw[:,1]
        self.sXYvert = raw[:,2]
        self.sYYvert = raw[:,3]
        #Get Horizontal Symmetry Plane Data
        raw = np.loadtxt(self.datapath + 'downPatch_sigmaxx_sigmaxy_sigmayy.xy')
        self.x = raw[:,0]
        self.sXXhorz = raw[:,1]
        self.sXYhorz = raw[:,2]
        self.sYYhorz = raw[:,3]
        #Get Diagonal Plane Data
        raw = np.loadtxt(self.datapath + 'diagPatch_sigmaxx_sigmaxy_sigmayy.xy')
        self.xdiag = raw[:,0]
        self.ydiag = raw[:,1]
        self.sXXdiag = raw[:,3]
        self.sXYdiag = raw[:,4]
        self.sYYdiag = raw[:,5]

def PlotXXVert(cases, name, show_plot):
    save_name = 'cs3.' + name + '.VertSym.sigXX.pdf'
    plt.figure()
    plt.title('X-Direction Normal Stress in Vertical Symmetry Plane', fontsize=txt_tit)
    plt.xlabel('$y\,[m]$', fontsize=txt_lbl)
    plt.ylabel('$(\sigma_{xx})_{x=0}\,[kPa]$', fontsize=txt_lbl)
    s = np.linspace(R, 2, N)
    sigma_analytic = VertSigmaXX(s)
    plt.plot(s, sigma_analytic, '-', label='Analyitic Solution', linewidth=line, markersize=mark)
    for i, case in enumerate(cases):
        #Case Name
        casename = str(case.L) + '_' + str(case.dx2) + '_case'
        label = 'FV Solution (L=' + str(case.L*2) + 'm, n$_{2}$=' + str(case.dx2) + ')'
        plt.plot(case.y, case.sXXvert, '-.', label=label, linewidth=line, markersize=mark)
    plt.legend(loc='best')
    plt.xlim(R,2)
    #Save Plot
    try:
        os.mkdir('Results')
    except Exception:
        pass
    plt.savefig('Results/' + save_name, bbox_inches='tight')
    if show_plot==1:
        plt.show()

def PlotYYVert(cases, name, show_plot):
    save_name = 'cs3.' + name + '.VertSym.sigYY.pdf'
    plt.figure()
    plt.title('Y-Direction Normal Stress in Vertical Symmetry Plane', fontsize=txt_tit)
    plt.xlabel('$y\,[m]$', fontsize=txt_lbl)
    plt.ylabel('$(\sigma_{yy})_{x=0}\,[kPa]$', fontsize=txt_lbl)
    s = np.linspace(R, 2, N)
    sigma_analytic = VertSigmaYY(s)
    plt.plot(s, sigma_analytic, '-', label='Analyitic Solution', linewidth=line, markersize=mark)
    for i, case in enumerate(cases):
        #Case Name
        casename = str(case.L) + '_' + str(case.dx2) + '_case'
        label = 'FV Solution (L=' + str(case.L*2) + 'm, n$_{2}$=' + str(case.dx2) + ')'
        plt.plot(case.y, case.sYYvert, '-.', label=label, linewidth=line, markersize=mark)
    plt.legend(loc='best')
    plt.xlim(R,2)
    #Save Plot
    try:
        os.mkdir('Results')
    except Exception:
        pass
    plt.savefig('Results/' + save_name, bbox_inches='tight')
    if show_plot==1:
        plt.show()

def PlotXXHorz(cases, name, show_plot):
    save_name = 'cs3.' + name + '.HorzSym.sigXX.pdf'
    plt.figure()
    plt.title('X-Direction Normal Stress in Horizontal Symmetry Plane', fontsize=txt_tit)
    plt.xlabel('$x\,[m]$', fontsize=txt_lbl)
    plt.ylabel('$(\sigma_{xx})_{y=0}\,[kPa]$', fontsize=txt_lbl)
    s = np.linspace(R, 2, N)
    sigma_analytic = HorzSigmaXX(s)
    plt.plot(s, sigma_analytic, '-', label='Analyitic Solution', linewidth=line, markersize=mark)
    for i, case in enumerate(cases):
        #Case Name
        casename = str(case.L) + '_' + str(case.dx2) + '_case'
        label = 'FV Solution (L=' + str(case.L*2) + 'm, n$_{2}$=' + str(case.dx2) + ')'
        plt.plot(case.x, case.sXXhorz, '-.', label=label, linewidth=line, markersize=mark)
    plt.legend(loc='best')
    plt.xlim(R,2)
    #Save Plot
    try:
        os.mkdir('Results')
    except Exception:
        pass
    plt.savefig('Results/' + save_name, bbox_inches='tight')
    if show_plot==1:
        plt.show()

def PlotYYHorz(cases, name, show_plot):
    save_name = 'cs3.' + name + '.HorzSym.sigYY.pdf'
    plt.figure()
    plt.title('Y-Direction Normal Stress in Horizontal Symmetry Plane', fontsize=txt_tit)
    plt.xlabel('$x\,[m]$', fontsize=txt_lbl)
    plt.ylabel('$(\sigma_{yy})_{y=0}\,[kPa]$', fontsize=txt_lbl)
    s = np.linspace(R, 2, N)
    sigma_analytic = HorzSigmaYY(s)
    plt.plot(s, sigma_analytic, '-', label='Analyitic Solution', linewidth=line, markersize=mark)
    for i, case in enumerate(cases):
        #Case Name
        casename = str(case.L) + '_' + str(case.dx2) + '_case'
        label = 'FV Solution (L=' + str(case.L*2) + 'm, n$_{2}$=' + str(case.dx2) + ')'
        plt.plot(case.x, case.sYYhorz, '-.', label=label, linewidth=line, markersize=mark)
    plt.legend(loc='best')
    plt.xlim(R,2)
    #Save Plot
    try:
        os.mkdir('Results')
    except Exception:
        pass
    plt.savefig('Results/' + save_name, bbox_inches='tight')
    if show_plot==1:
        plt.show()

def PlotXXDiag(cases, name, show_plot):
    save_name = 'cs3.' + name + '.Diag.sigXX.pdf'
    plt.figure()
    plt.title('X-Direction Normal Stress in Diagonal Plane', fontsize=txt_tit)
    plt.xlabel('$x,y\,[m]$', fontsize=txt_lbl)
    plt.ylabel('$(\sigma_{xx})_{x=y}\,[kPa]$', fontsize=txt_lbl)
    s = np.linspace(0.707107, 2, N)
    sigma_analytic = DiagSigmaXX(s)
    plt.plot(s, sigma_analytic, '-', label='Analyitic Solution', linewidth=line, markersize=mark)
    for i, case in enumerate(cases):
        #Case Name
        casename = str(case.L) + '_' + str(case.dx2) + '_case'
        label = 'FV Solution (L=' + str(case.L*2) + 'm, n$_{2}$=' + str(case.dx2) + ')'
        plt.plot(case.xdiag, case.sXXdiag, '-.', label=label, linewidth=line, markersize=mark)
    plt.legend(loc='best')
    plt.xlim(R,2)
    #Save Plot
    try:
        os.mkdir('Results')
    except Exception:
        pass
    plt.savefig('Results/' + save_name, bbox_inches='tight')
    if show_plot==1:
        plt.show()

def PlotYYDiag(cases, name, show_plot):
    save_name = 'cs3.' + name + '.Diag.sigYY.pdf'
    plt.figure()
    plt.title('Y-Direction Normal Stress in Diagonal Plane', fontsize=txt_tit)
    plt.xlabel('$x,y\,[m]$', fontsize=txt_lbl)
    plt.ylabel('$(\sigma_{yy})_{x=y}\,[kPa]$', fontsize=txt_lbl)
    s = np.linspace(0.707107, 2, N)
    sigma_analytic = DiagSigmaYY(s)
    plt.plot(s, sigma_analytic, '-', label='Analyitic Solution', linewidth=line, markersize=mark)
    for i, case in enumerate(cases):
        #Case Name
        casename = str(case.L) + '_' + str(case.dx2) + '_case'
        label = 'FV Solution (L=' + str(case.L*2) + 'm, n$_{2}$=' + str(case.dx2) + ')'
        plt.plot(case.xdiag, case.sYYdiag, '-.', label=label, linewidth=line, markersize=mark)
    plt.legend(loc='best')
    plt.xlim(R,2)
    #Save Plot
    try:
        os.mkdir('Results')
    except Exception:
        pass
    plt.savefig('Results/' + save_name, bbox_inches='tight')
    if show_plot==1:
        plt.show()

def RMSerror(numeric, analytic):
    n = len(numeric)
    rms = 0
    for i in range(0,n):
        rms+= ((numeric[i] - analytic[i]) ** 2.)

    rms = (rms / n) ** 0.5
    return rms

def CaseRMS(cases, name, show_plot):
    save_name = 'cs3.' + name + '.RMS.pdf'
    RMSxxVert = np.zeros(len(cases))
    RMSyyVert = np.zeros(len(cases))
    RMSyyHorz = np.zeros(len(cases))
    RMSxxDiag = np.zeros(len(cases))
    RMSxxVertMax = 0
    RMSyyVertMax = 0
    for i, case in enumerate(cases):
        if np.amax(case.sXXvert) > RMSxxVertMax:
            RMSxxVertMax = np.amax(case.sXXvert)
        if np.amax(case.sYYvert) > RMSyyVertMax:
            RMSyyVertMax = np.amax(case.sYYvert)

    for i, case in enumerate(cases):
        RMSxxVert[i] = RMSerror(case.sXXvert, VertSigmaXX(case.y))/sigma
        RMSyyVert[i] = RMSerror(case.sYYvert, VertSigmaYY(case.y))/sigma
        RMSyyHorz[i] = RMSerror(case.sYYhorz, HorzSigmaYY(case.x))/sigma
        RMSxxDiag[i] = RMSerror(case.sXXdiag, HorzSigmaYY(case.xdiag))/sigma

    print(RMSxxVert, RMSyyHorz, RMSxxDiag)
    plt.figure()
    plt.title('Normalized RMS Error', fontsize=txt_tit)
    plt.xlabel('$n_{2}$', fontsize=txt_lbl)
    plt.ylabel('$RMS/\sigma [n.d]$', fontsize=txt_lbl)
    plt.plot([c.dx2 for c in cases], RMSxxVert, '-.', label='$(\sigma_{xx})_{x=0}$')
    # plt.plot([c.dx2 for c in cases], RMSyyVert, '-.', label='$(\sigma_{yy})_{x=0}$')
    plt.legend(loc='best')
    #Save Plot
    try:
        os.mkdir('Results')
    except Exception:
        pass
    plt.savefig('Results/' + save_name, bbox_inches='tight')
    if show_plot==1:
        plt.show()

def main(lengths, spacing2, name):
    global R, sigma, N
    R = 0.5             #perforation radius
    sigma = 10000       #applied stress in x-dir
    N = 101
    #Create data classes
    n_cases = len(lengths)
    cases = np.zeros(n_cases, dtype=object)
    for i, (L, dx2) in enumerate(zip(lengths, spacing2)):
        cases[i] = CaseData(L, dx2)
    PlotXXVert(cases, name, 0)
    PlotYYVert(cases, name, 0)
    PlotXXHorz(cases, name, 0)
    PlotYYHorz(cases, name, 0)
    PlotXXDiag(cases, name, 0)
    PlotYYDiag(cases, name, 0)
    CaseRMS(cases, name, 0)
    # PlotSigma(cases, 'vert', 'xx', 'Test', 1)


if __name__ == "__main__":
    
    # name = 'MeshSens'
    # lengths =  [2, 2,  2,  2 , 2,  2]
    # spacing1 = [5, 10, 20, 40, 80, 200]
    # spacing2 = 2*spacing1
    # main(lengths, spacing2, name)

    # name = 'MeshSens'
    # lengths =  [2,  2,  2, 2, 2]
    # spacing2 = [5, 10, 40, 200, 500]

    # main(lengths, spacing2, name)

    name = 'Length'
    lengths =  [2,  50, 100]
    spacing2 = [10, 500, 1000]
    main(lengths, spacing2, name)