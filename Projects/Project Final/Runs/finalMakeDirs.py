"""
FINAL - CHANNEL FLOW WITH CAVITY - DIRECTORY SETUP
MAE 219
LOGAN HALSTROM
15 DECEMBER 2014

DESCRIPTION:  This code automates OpenFoam directory setup.
"""

import numpy as np
import matplotlib.pyplot as plt 
import subprocess
import os

def cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()

def MakeMeshFile(case, AR, nx_cav, ny_cav, nx_cha, ny_cha, sx, sy):
    """Simple case: uniform(ish) mesh in cavity, graded mesh in channel. 
    sx, sy expansion ratios of channel 
    """
    
    towrite = '''
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.3.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      blockMeshDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    convertToMeters 1;
    vertices
    (
        (''' + str(-L/2) + ''' 0 0)
        (''' + str(-L/2) + ''' ''' + str(-AR*L) + ''' 0)
        (''' + str(L/2) + ''' ''' + str(-AR*L) + ''' 0)
        (''' + str(L/2) + ''' 0 0)
        (''' + str(chan_len/2) + ''' 0 0)
        (''' + str(chan_len/2) + ''' ''' + str(L) + ''' 0)
        (''' + str(L/2) + ''' ''' + str(L) + ''' 0)
        (''' + str(-L/2) + ''' ''' + str(L) + ''' 0)
        (''' + str(-chan_len/2) + ''' ''' + str(L) + ''' 0)
        (''' + str(-chan_len/2) + ''' 0 0)
        (''' + str(-L/2) + ''' 0 0.1)
        (''' + str(-L/2) + ''' ''' + str(-AR*L) + ''' 0.1)
        (''' + str(L/2) + ''' ''' + str(-AR*L) + ''' 0.1)
        (''' + str(L/2) + ''' 0 0.1)
        (''' + str(chan_len/2) + ''' 0 0.1)
        (''' + str(chan_len/2) + ''' ''' + str(L) + ''' 0.1)
        (''' + str(L/2) + ''' ''' + str(L) + ''' 0.1)
        (''' + str(-L/2) + ''' ''' + str(L) + ''' 0.1)
        (''' + str(-chan_len/2) + ''' ''' + str(L) + ''' 0.1)
        (''' + str(-chan_len/2) + ''' 0 0.1)

    );
    blocks
    (
        hex (1 2 3 0 11 12 13 10) (''' + str(int(nx_cav)) + ' ' + str(int(ny_cav)) + ''' 1) simpleGrading (1 1 1)
        hex (3 4 5 6 13 14 15 16) (''' + str(int(nx_cha)) + ' ' + str(int(ny_cha)) + ''' 1) simpleGrading (''' + str(sx) + ' ' + str(sy) + ''' 1)
        hex (0 3 6 7 10 13 16 17) (''' + str(int(nx_cav)) + ' ' + str(int(ny_cha)) + ''' 1) simpleGrading (1 ''' + str(sy) + ''' 1)
        hex (9 0 7 8 19 10 17 18) (''' + str(int(nx_cha)) + ' ' + str(int(ny_cha)) + ''' 1) simpleGrading (''' + str(1./sx) + ' ' + str(sy) + ''' 1)
    );
    edges
    (
    );
    boundary
    (
        inlet
        {
            type patch;
            faces
            (
                (9 8 18 19)
            );
        }
        outlet
        {
            type patch;
            faces
            (
                (4 14 15 5)
            );
        }
        movingWalls
        {
            type wall;
            faces
            (
                (5 15 16 6)
                (6 16 17 7)
                (7 17 18 8)
            );
        }
        fixedWalls
        {
            type wall;
            faces
            (
                (3 13 14 4)
                (2 12 13 3)
                (1 11 12 2)
                (1 0 10 11)
                (9 19 10 0)
            );
        }
        frontAndBack
        {
            type empty;
            faces
            (
                (3 4 5 6)
                (0 3 6 7)
                (9 0 7 8)
                (1 2 3 0)
                (13 16 15 14)
                (10 17 16 13)
                (19 18 17 10)
                (11 10 13 12)
            );
        }
    );
    mergePatchPairs
    (
    );
    // ************************************************************************* //
    '''

    #Rewrite blockMeshDict to create correct mesh
    path = case +  '/constant/polyMesh/blockMeshDict'
    with open(path, 'w') as writefile:
        writefile.write(towrite)

def MakeTransportProps(case, Re, l):
    """Write transportProperties.  Calculate appropriate kinematic viscoity based on
    Reynold's number, characteristic length, and flow velocity.
    """

    nu = l * u / Re    #kinematic viscosity [m^2/s]    
    towrite1 = '''
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.3.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      transportProperties;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    '''
    
    if solver == 'pisoFoam':
        towrite2 = '''
        transportModel Newtonian;'''
    else: 
        towrite2 = ''

    towrite3 = '''

    nu              nu [ 0 2 -1 0 0 0 0 ] ''' + str(nu) + ''';


    // ************************************************************************* //
    '''
    towrite = towrite1 + towrite2 + towrite3
    #Rewrite blockMeshDict to create correct mesh
    path = case +  '/constant/transportProperties'
    with open(path, 'w') as writefile:
        writefile.write(towrite)

def MakePisoFiles(case, epsilon, k, nut):
    """Write files for pisoFoam solver: epsilon, k, nut
    """

    #WRITE EPSILON FILE
    towrite = '''
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.3.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        location    "0";
        object      epsilon;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 2 -3 0 0 0 0];

    internalField   uniform ''' + str(epsilon) + ''';

    boundaryField
    {
        inlet
        {
            type            zeroGradient;
        }
        outlet
        {
            type            zeroGradient;
        }
        movingWalls
        {
            type            epsilonWallFunction;
            value           uniform ''' + str(epsilon) + ''';
        }
        fixedWalls
        {
            type            epsilonWallFunction;
            value           uniform ''' + str(epsilon) + ''';
        }
        frontAndBack
        {
            type            empty;
        }
    }


    // ************************************************************************* //
    '''
    #Rewrite 
    path = case +  '/0/epsilon'
    with open(path, 'w') as writefile:
        writefile.write(towrite)

    #WRITE k FILE
    towrite = '''
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.3.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        location    "0";
        object      k;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 2 -2 0 0 0 0];

    internalField   uniform ''' + str(k) + ''';

    boundaryField
    {
        inlet
        {
            type            zeroGradient;
        }
        outlet
        {
            type            zeroGradient;
        }
        movingWalls
        {
            type            kqRWallFunction;
            value           uniform ''' + str(k) + ''';
        }
        fixedWalls
        {
            type            kqRWallFunction;
            value           uniform ''' + str(k) + ''';
        }
        frontAndBack
        {
            type            empty;
        }
    }
    '''
    #Rewrite 
    path = case +  '/0/k'
    with open(path, 'w') as writefile:
        writefile.write(towrite)

    #WRITE nut FILE
    towrite = '''
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.3.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        location    "0";
        object      nut;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 2 -1 0 0 0 0];

    internalField   uniform ''' + str(nut) + ''';

    boundaryField
    {
        inlet
        {
            type            zeroGradient;
        }
        outlet
        {
            type            zeroGradient;
        }
        movingWalls
        {
            type            nutkWallFunction;
            value           uniform ''' + str(nut) + ''';
        }
        fixedWalls
        {
            type            nutkWallFunction;
            value           uniform ''' + str(nut) + ''';
        }
        frontAndBack
        {
            type            empty;
        }
    }
    '''
    #Rewrite 
    path = case +  '/0/nut'
    with open(path, 'w') as writefile:
        writefile.write(towrite)

def MakeControlDict(case, solver, nx, t_start, t_end, nsave):
    """Write controlDict.  Input solver name as string, number of points in cavity
    to determine time step based of CFL condition, end time for simulation, and 
    number of intervals to save
    """

    dx = L / nx
    dt = CFL * dx / u
    dsave = int(((t_end - t_start) / nsave) / dt)       #save every dsave intervals 
    if dsave < 1:
        dsave = 1   
    towrite = '''
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.3.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      controlDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    application     ''' + solver + ''';

    startFrom       startTime;

    startTime       ''' + str(t_start) + ''';

    stopAt          endTime;

    endTime         ''' + str(t_end) + ''';

    deltaT          ''' + str(dt) + ''';

    writeControl    timeStep;

    writeInterval   ''' + str(dsave) + ''';

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;


    // ************************************************************************* //
    '''

    #Rewrite blockMeshDict to create correct mesh
    path = case +  '/system/controlDict'
    with open(path, 'w') as writefile:
        writefile.write(towrite)

def MakeSampleDict(case, xupstream, AR):
    """Write sampleDict.  Sample upstream velocity profile to guage convergence.
    Sample centerline velocity profile.
    xupstream --> location to sample upstream of centerline
    """

    towrite = '''
    
    /*--------------------------------*- C++ -*----------------------------------*    | =========                 |                                                 |
    | \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \    /   O peration     | Version:  2.3.0                                 |
    |   \  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      sampleDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    interpolationScheme cellPoint;

    setFormat       raw;

    sets
    (
        centerLine
        {
            type    uniform;
            axis    y;
            start   ( 0 ''' + str(-AR*L) + ''' 0.05 );
            end     ( 0 1.0 0.05 );
            nPoints 100;
        }
        freeStream
        {
            type    uniform;
            axis    y;
            start   ( ''' + str(xupstream) + ''' 0.0 0.05 );
            end     ( ''' + str(xupstream) + ''' 1.0 0.05 );
            nPoints 100;
        }
    );
    
    fields
    (
        Ux
        Uy
    );


    // ************************************************************************* //
    '''

    #Rewrite blockMeshDict to create correct mesh
    path = case +  '/system/sampleDict'
    with open(path, 'w') as writefile:
        writefile.write(towrite)

def MakeCaseDir(case):
    """ Given arrays of domain lengths and cell spacings to run,
    copy base directory to new folders specific for each case.
    case --> path of case directory
    """
    #Copy base directory to case name if it does not already exist
    if not os.path.exists(case):
        command = 'cp -r baseCase/ ' + case + '/; '
        cmd(command)

#DOMAIN PARAMETERS
Re = 100.             #Reynold's number
AR = 2            #Aspect Ratio of cavity (h/L)
L = 1.              #Length of cavity
chan_len = L*10     #Length of channel (approx. infinity)
H = L               #Height of channel
h = AR * L          #Height of cavity

u = 1.              #freestream velocity
CFL = 1.            #CFL number
t_start = 60
t_end = 60.1
nsave = 10              #number of save intervals in solution
solver = 'icoFoam'     #solver
epsilon = 0.000765
k = 0.00325
nut = 0

#GRID PARAMETERS
# #TRIAL RUNS**************************************
# #Cavity cell size
# nx_cav = 50
# #Channel cell size
# nx_cha = 60
# ny_cha = 20
# sx, sy = 7, 4 #expansion ratios for channel
#ULTRAFINE
nx_cav = 1000
nx_cha = 60
ny_cha = 20
sx, sy = 7, 4 #expansion ratios for channel
#OTHER
# #Cavity cell size
# nx_cav = 100
# #Channel cell size
# nx_cha = 70
# ny_cha = 25
# sx, sy = 8, 6 #expansion ratios for channel
# #FINER
# nx_cav = 60
# nx_cha = 70
# ny_cha = 25
# sx, sy = 5, 5

ny_cav = nx_cav * h / L     #same spacing x and y for cavity
dx_cav = L / nx_cav
dy_cav = AR * L / ny_cav
# dy_cav = dx_cav
# nx_cav = dx_cav * L
# ny_cav = dy_cav * h

#CASE NAME
tag = ''      #tag for case directory name (lead with underscore)
case = 'ar%s_Re%1.0f_dx%1.3f'%(AR, Re, dx_cav) + tag
#MAKE CASE DIRECTORY
MakeCaseDir(case)
#REWRITE blockMeshDict
MakeMeshFile(case, AR, nx_cav, ny_cav, nx_cha, ny_cha, sx, sy)
#REWRITE transportProperties
MakeTransportProps(case, Re, L)
if solver == 'pisoFoam':
    MakePisoFiles(case, epsilon, k, nut)
#REWRITE controlDict
MakeControlDict(case, solver, nx_cav, t_start, t_end, nsave)
MakeSampleDict(case, -2.5, AR)

#INFORMATION OUTPUT
dt = CFL * dx_cav / u
print('Case: AR=%s, Re=%s'%(AR, Re))
print('Cavity Grid: nx=%s, ny=%s, dx=%s, dy=%s, dt=%s'%(nx_cav, ny_cav, dx_cav, dy_cav, dt))
print('Steps to Run: ' + str(int((t_end - t_start) / dt)))



