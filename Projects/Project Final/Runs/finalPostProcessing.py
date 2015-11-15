"""
FINAL - CHANNEL FLOW WITH CAVITY - POST PROCESSING
MAE 219
LOGAN HALSTROM
16 DECEMBER 2014

DESCRIPTION:  This code automates OpenFoam directory setup.
"""

import numpy as np
import matplotlib.pyplot as plt 
import subprocess
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

#Line Styles
mark = 8
line = 1.5
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
filetype = '.pdf'

def MakeOutputDir(Dir):
    """make results output directory if it does not already exist.  Input directory name"""
    try:
        os.mkdir(Dir)
    except Exception:
        pass

def ReadVelComps(case, loc, tend):
    """
    case --> case directory path
    loc --> 'centerLine' or 'freeStream'
    """
    tplot = str(tend)     #time to plot (folder name)     
    path = case + '/postProcessing/sets/' + tplot + '/' + loc + '_Ux_Uy.xy'
    #DATA FROM VERTICAL LINE AT CENTERLINE OR OFFSET AHEAD IN THE FREESTREAM
    raw = np.loadtxt(path)
    y = raw[:,0]
    u = raw[:,1]
    v = raw[:,2]

    return y, u, v

def PlotFreestreamVel(case, tend, dt, showplot):
    """Plot convergence of velocity profile in freestream
    """
    plt.figure()
    plt.title('Freestream X-Velocity Convergence\nAR=%s, Re=%s, x=-2.5'%(AR, Re), fontdict=font_tit)
    plt.xlabel(r'$u$', fontdict=font_lbl)
    plt.ylabel(r'$y$', fontdict=font_lbl)
    for i in range(0, tend/dt+1):
        tplot = str(i*dt)     #time to plot (folder name)     
        path = case + '/postProcessing/sets/' + tplot + '/' + 'freeStream' + '_Ux_Uy.xy'
        #DATA FROM VERTICAL LINE AT CENTERLINE OR OFFSET AHEAD IN THE FREESTREAM
        raw = np.loadtxt(path)
        y = raw[:,0]
        u = raw[:,1]
        v = raw[:,2]
        plt.plot(u, y, '-', label='t=' + tplot + 's')

    plt.legend(loc='best', fancybox=True, framealpha=0.5)
    save_name = case + '_FreestreamConverge' + filetype
    plt.savefig(ouputdir + '/' + save_name, bbox_inches='tight')
    if showplot == 1:
        plt.show()



#SINGLE CASE PARAMETERS
AR = 0.5
Re = 100
dx = 0.020
case = 'ar%s_Re%1.0f_dx%1.3f'%(AR, Re, dx)
tplot = 60

#RESULTS OUTPUT DIRECTORY
ouputdir = 'Results'
MakeOutputDir(ouputdir)

#PLOT FREESTREAM VELOCITY PROFILE AS IT CONVERGES
incr = 4
PlotFreestreamVel(case, tplot, incr, 0)

#PLOT X-VELOCITY ALONG CENTERLINE
# location = 'freeStream'
location = 'centerLine'
y, u, v = ReadVelComps(case, location, tplot)
plt.figure()
if location == 'centerLine':
    loc = 'Centerline'
else:
    loc = 'Freestream'

plt.title(loc + ' X-Velocity Profile\nAR=%s, Re=%s, t=%ss'%(AR, Re, tplot), fontdict=font_tit)
plt.xlabel(r'$u\,[m/s]$', fontdict=font_lbl)
plt.ylabel(r'$y\,[m]$', fontdict=font_lbl)
line = 1.5
vertlinex = np.zeros(len(y))
plt.plot(vertlinex, y, 'g', linewidth=line)
plt.fill_betweenx(y, vertlinex, u, facecolor='green', alpha=0.2)
wd, ln = 0.03, 0.03
for i in range(0, len(y), 4):
    if abs(u[i]) < ln:
        plt.plot([0, u[i]], [y[i], y[i]], 'g', linewidth=line)
    else:
        plt.arrow(0, y[i], u[i]-ln, 0, head_width=wd, head_length=ln,
            fc='g', ec='g', linewidth=line)
plt.plot(u, y, 'g', linewidth=line)
plt.axis([min(u), max(u), min(y), max(y)])
save_name = case + '_' + loc +'Vel' + filetype
plt.savefig(ouputdir + '/' + save_name, bbox_inches='tight')
plt.show()
