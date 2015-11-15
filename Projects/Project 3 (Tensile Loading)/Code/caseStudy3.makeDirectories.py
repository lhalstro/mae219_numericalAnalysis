"""
CASE STUDY 3 - PERFORATED PLATE IN TENSION
MAE 219
LOGAN HALSTROM
15 NOVEMBER 2014

DESCRIPTION:  This code automates OpenFoam simulations of the test case
for given case conditions and then produces plots of stresses compared
to the analytic result.
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

def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()

def MakeMeshFile(L, dx1, dx2):
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
        (0.5 0 0)
        (1 0 0)
        (''' + str(L) + ''' 0 0)
        (''' + str(L) + ''' 0.707107 0)
        (0.707107 0.707107 0)
        (0.353553 0.353553 0)
        (''' + str(L) + ''' 2 0)
        (0.707107 2 0)
        (0 2 0)
        (0 1 0)
        (0 0.5 0)
        (0.5 0 0.5)
        (1 0 0.5)
        (''' + str(L) + ''' 0 0.5)
        (''' + str(L) + ''' 0.707107 0.5)
        (0.707107 0.707107 0.5)
        (0.353553 0.353553 0.5)
        (''' + str(L) + ''' 2 0.5)
        (0.707107 2 0.5)
        (0 2 0.5)
        (0 1 0.5)
        (0 0.5 0.5)
    );
    blocks
    (
        hex (5 4 9 10 16 15 20 21) (''' + str(dx1) + ' ' + str(dx1) + ''' 1) simpleGrading (1 1 1)
        hex (0 1 4 5 11 12 15 16) (''' + str(dx1) + ' ' + str(dx1) + ''' 1) simpleGrading (1 1 1)
        hex (1 2 3 4 12 13 14 15) (''' + str(dx2) + ' ' + str(dx1) + ''' 1) simpleGrading (1 1 1)
        hex (4 3 6 7 15 14 17 18) (''' + str(dx2) + ' ' + str(dx1 *2) + ''' 1) simpleGrading (1 1 1)
        hex (9 4 7 8 20 15 18 19) (''' + str(dx1) + ' ' + str(dx1 * 2) + ''' 1) simpleGrading (1 1 1)
    );
    edges
    (
        arc 0 5 (0.469846 0.17101 0)
        arc 5 10 (0.17101 0.469846 0)
        arc 1 4 (0.939693 0.34202 0)
        arc 4 9 (0.34202 0.939693 0)
        arc 11 16 (0.469846 0.17101 0.5)
        arc 16 21 (0.17101 0.469846 0.5)
        arc 12 15 (0.939693 0.34202 0.5)
        arc 15 20 (0.34202 0.939693 0.5)
    );
    boundary
    (
        left
        {
            type symmetryPlane;
            faces
            (
                (8 9 20 19)
                (9 10 21 20)
            );
        }
        right
        {
            type patch;
            faces
            (
                (2 3 14 13)
                (3 6 17 14)
            );
        }
        down
        {
            type symmetryPlane;
            faces
            (
                (0 1 12 11)
                (1 2 13 12)
            );
        }
        up
        {
            type patch;
            faces
            (
                (7 8 19 18)
                (6 7 18 17)
            );
        }
        hole
        {
            type patch;
            faces
            (
                (10 5 16 21)
                (5 0 11 16)
            );
        }
        frontAndBack
        {
            type empty;
            faces
            (
                (10 9 4 5)
                (5 4 1 0)
                (1 4 3 2)
                (4 7 6 3)
                (4 9 8 7)
                (21 16 15 20)
                (16 11 12 15)
                (12 13 14 15)
                (15 14 17 18)
                (15 18 19 20)
            );
        }
    );
    mergePatchPairs
    (
    );
    // ************************************************************************* //
    '''

    return towrite


def MakeSampleFile(L):
    """write simpleDict file to calculate normal stresses
    at planes of symmetry"""
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
        object      sampleDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    interpolationScheme cellPoint;

    setFormat       raw;

    sets
    (
        leftPatch
        {
            type    uniform;
            axis    y;
            start   ( 0 0.5 0.25 );
            end     ( 0 2 0.25 );
            nPoints 100;
        }
        downPatch
        {
            type    uniform;
            axis    x;
            start   ( 0.5 0 0.25 );
            end     ( ''' + str(L) + ''' 0 0.25 );
            nPoints ''' + str(int(100*L/2)) + ''';
        }
        diagPatch
        {
            type    uniform;
            axis    xyz;
            start   ( 0.707107 0.707107 0.25 );
            end     ( ''' + str(L) + ''' 2 0.25 );
            nPoints ''' + str(int(100*L/2)) + ''';
        }
    );

    fields
    (
        sigmaxx
        sigmayy
        sigmaxy
    );


    // ************************************************************************* //
    '''

    return towrite

def MakeCaseDirs(length, spacing1, spacing2):
    """ Given arrays of domain lengths and cell spacings to run,
    copy base directory to new folders specific for each case.
    length --> Plate lengths
    spacing1 --> number of cells in x for block 1
    spacing2 --> number of cells in x for block 2
    """
    for  L, dx1, dx2 in zip(length, spacing1, spacing2):
        #Case Name
        case = str(L) + '_' + str(dx2) + '_case' 
        #Copy base directory to case name if it does not already exist
        if not os.path.exists(case):
            command = 'cp -r plateHole/ ' + case + '/; '
            subprocess_cmd(command)

        #Rewrite blockMeshDict to create correct mesh
        path = case +  '/constant/polyMesh/blockMeshDict'
        with open(path, 'w') as writefile:
            writefile.write(MakeMeshFile(L, dx1, dx2))

        #Rewrite sampleDict to calculate correct symmetry plane stress
        path = case +  '/system/sampleDict'
        with open(path, 'w') as writefile:
            writefile.write(MakeSampleFile(L))

 #        print ('Mesh file generated')

 #      #Run Cases
 #          #if the t=100 save doesn't exist, run
 #      if not os.path.exists(case + '/100/'):
 #            print(case + ' running')
 #            # command = "mount $HOME/OpenFOAM OpenFOAM.sparsebundle; "
 #            # command = "of230; "
 #            # command = "source $HOME/OpenFOAM/OpenFOAM-2.3.0/etc/bashrc WM_NCOMPPROCS=6 WM_MPLIB=SYSTEMOPENMPI; "
 #            #Enter case directory
 #            command = "hdiutil attach -quiet -mountpoint $HOME/OpenFOAM OpenFOAM.sparsebundle; "
 #            command += "sleep 1; "
 #            command += "cd " + case + '; '
 #            #Generate Mesh
 #            command += "blockMesh; "
 #            #Run Case
 #            command += "solidDisplacementFoam > log; "
 #            #Calculate Stresses
 #            command += "foamCalc components sigma; "
 #            #Find stress in line
 #            command += "sample"
 #            subprocess_cmd(command)
 #        print(case + ' complete.')

    # print('Simulations complete.')


def main(length, spacing1, spacing2):
    """ Run simulation for given plate lengths and cell spacings
    length --> Plate lengths
    spacing1 --> number of cells in x for block 1
    spacing2 --> number of cells in x for block 2
    """
    MakeCaseDirs(length, spacing1, spacing2)


if __name__ == "__main__":
    # # Test Case
    # length =   [2,  2]
    # spacing1 = [10, 20]
    # spacing2 = [20, 40]
    # main(length, spacing1, spacing2)

    #Mesh Sensitivity Study
    lengths =  [2, 2,  2,  2 , 2,  2,   2,   2]
    spacing1 = [5, 10, 20, 40, 80, 200, 500]
    spacing2 = 2*spacing1
    main(lengths, spacing1, spacing2)

    # #Length Sensitivity
    # lengths =  [2,  10,    50,    100]
    # spacing1 = [10, 10, 10, 10]
    # spacing2 = 2*[10, 10*10, 10*50, 10*100]
    # main(lengths, spacing1, spacing2)