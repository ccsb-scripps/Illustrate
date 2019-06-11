[![Build Status](https://api.travis-ci.org/ccsb-scripps/Illustrate.svg?branch=master)](https://api.travis-ci.org/ccsb-scripps/Illustrate.svg?branch=master)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
**I L L U S T R A T E**
Biomolecular Illustration
copyright 2019 David S Goodsell

Licensed under the Apache License, Version 2.0 (the "License") you may not use this file except in compliance with the License. You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License

This work was supported in part by Damon Runyon-Walter Winchell Cancer Research Fund Fellowship DRG 972, the US National Institutes of Health R01-GM120604 and the kind support of the RCSB Protein Data Bank.
DS Goodsell & AJ Olson (1992) "Molecular Illustration in Black and White" JMolGraphics 10, 235-240.
April 2019 Simplified and released with only non-photorealistic rendering

**COMPILING AND RUNNING**
Compile the fortran code:

    gfortran illustrate.f -o illustrate

Run the program:

    illustrate < command_file

Command file is read from unit 5 (standard in), and a bunch of diagnostic stuff is written to unit 6 (standard out)

**COMMAND FILE FORMAT**
The command file has command cards (read, center, world, calculate, etc), followed by parameter cards needed for each command. Please issue command cards in this order:
1. `read` — reads coordinates and selection/rendering parameters
2. `center, translate, xrot, yrot, zrot, scale` commands, in any order
3. `world` — defines rendering parameter
4. `illustrate` — defines illustration parameters
5. `calculate` — renders the image and writes ppm file

This program has many idiosyncracies and almost no error checking (beware…postdoc code from the 90s!)

Selection/rendering cards in the “read” command are read sequentially until one without “ATOM” or “HETATM” is found. Atoms are compared to cards in order, and if a match is found, the atom is assigned those parameters. This process provides a lot of flexibility with very few cards, if you’re clever about the order of cards and the use of wildcards. Any number of rotation cards may be added, and they are concatenated when added. This means they are effectively applied last to first, so if you’re progressively refining an orientation, add new rotations to the top of the list. Rotations are applied first, then centering, and finally translation. This ensures that the molecule is always centered in the view. Use the translation if you want it offset. Origin at upper left, +x down, +y left to right, +z towards viewer, molecules clipped at z=0

**ABOUT THE OUTLINES**
Outlines are created by calculating the local derivative of the z-value using a variety of kernels, then providing a range of values of these derivatives that will be given shades of gray and black. Outlines may be drawn based on the contours and outer shape of the molecule, between subunits, and based on the differences in the residue number in chains.

**SAMPLE COMMAND FILE** for PDB entry 2hhb (Hemoglobin)

    read                                        #READ command
    2hhb.pdb                                         #PDB format coordinate file
    HETATM-----HOH-- 0,9999, 0.5,0.5,0.5, 0.0        #selection/rendering cards
    ATOM  -H-------- 0,9999, 0.5,0.5,0.5, 0.0          # omits hydrogens
    ATOM  H--------- 0,9999, 0.5,0.5,0.5, 0.0          # omits hydrogens with long atom names
    ATOM  -C-------A 0,9999, 1.0,0.6,0.6, 1.6          # draws carbon atoms in chain A, pink
    ATOM  -S-------A 0,9999, 1.0,0.5,0.5, 1.8
    ATOM  ---------A 0,9999, 1.0,0.5,0.5, 1.5          # draws the rest of chain A atoms
    ATOM  -C-------C 0,9999, 1.0,0.6,0.6, 1.6          # similarly for chain C
    ATOM  -S-------C 0,9999, 1.0,0.5,0.5, 1.8
    ATOM  ---------C 0,9999, 1.0,0.5,0.5, 1.5
    ATOM  -C-------- 0,9999, 1.0,0.8,0.6, 1.6          # draws remaining protein carbons, light orange
    ATOM  -S-------- 0,9999, 1.0,0.7,0.5, 1.8
    ATOM  ---------- 0,9999, 1.0,0.7,0.5, 1.5
    HETATMFE---HEM-- 0,9999, 1.0,0.8,0.0, 1.8          # draws heme iron, yellow
    HETATM-C---HEM-- 0,9999, 1.0,0.3,0.3, 1.6
    HETATM-----HEM-- 0,9999, 1.0,0.1,0.1, 1.5
    END                                                #end of READ command
    center                                       #CENTER command
    auto                                             #use autocentering
    trans                                        #TRANSLATION command
    0.,0.,0.                                         #x,y,z for translation
    scale                                        #SCALE command
    12.0                                             #scale value (pixels/A)
    zrot                                         #ROTATION command
    90.                                              #rotation angle (deg)
    wor                                          #WORLD command (rendering parameters)
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0                  #rgbback, rgbfog, fogparams
    1,0.0023,2.0,1.0,0.2                             #soft shadow parameters
    -30,-30                                          #image size with autosizing
    illustrate                                   #ILLUSTRATION command
    3.0,10.0,4,0.0,5.0                               #contour outlines
    3.0,10.0                                         #subunit outlines
    3.0,8.0,6000.                                    #residue outlines
    calculate                                    #CALCULATE command
    2hhb.pnm                                         #image file name in ppm 
 


**COMMAND CARDS**

READ command

    PDB filename

    Selection/rendering cards

         record name (A6) “ATOM  “ or “HETATM”, matched with columns 1-6 of PDB file

         atom descriptor (A10) matched with columns 13-22 of PDB file, “-” is a wildcard

         residue range low, high (integer)

         color r,g,b (3 real)  0.0-1.0

         radius (real) angstrom, value of 0.0 will omit the atom

--------------------------------------------------------------------
CENTER command

    centering type (a3)
        “aut”  (auto) will autocenter with no clipping

        “cen” (center) will center on max/min coordinates in xyz

--------------------------------------------------------------------
TRANSLATE command

    translation vector x,y,z (3 real) translation in Angstroms

--------------------------------------------------------------------
XROT, YROT, ZROT commands

    rotation angle (real) rotation angle in degrees

--------------------------------------------------------------------
SCALE command

    scale value (real) pixels/Angstrom

--------------------------------------------------------------------
WORLD command

Background and fog parameters

    rback r,g,b (3 real) color of background (0.0-1.0)

    rfog r,g,b (3 real) color of fog (0.0-1.0)

    pfogh, pfogl (2 real) fractional transparency of fog at front and back of molecule

                                       (0.0-1.0, 1.0=no fog)

Soft shadow parameters

    icone (integer) 0=no shadows, 1=shadows

    pcone (real) fractional shadowing around each atom

                             larger=darker (0.0-1.0, typically 0.0023)

    coneangle (real) angle of shadowing around each atom

                            larger=tighter region (0.0-1.0, typically 2.0)

    rcone (real) shadowing only applied if z-difference greater than this value

                           (Angstroms, typically 1.0)

    pshadowmax (real) maximal shadowing amount

                            smaller=darker (0.0-1.0, typically 0.7)

--------------------------------------------------------------------
ILLUSTRATE command

Contour outline parameters

    l_low, l_high (2 real) thresholds for gray to black

                 typically values from about 3.0-20.0, best values for typical atomic illustrations: 3.0, 10.0

    ikernel (integer) kernel for derivative calculation (1,2,3,4 smoothest=4)

    l_diff_min, l_diff_max (2 integer) range of z-difference used for derivative (Angstroms)

                    0.0,1.0 gives outlines around every atom

                    0.0,1000.0 gives only outline around molecule

                    0.0,5.0 is typical

Subunit outline parameters

    l_low, l_high (2 real) thresholds for gray to black (typically ~ 3.0-20.0)

Residue outline parameters

    l_low, l_high (2 real) thresholds for gray to black (typically ~ 3.0-20.0)

    resdiff (real) difference in residue numbers to draw outlines

--------------------------------------------------------------------
CALCULATE command

    Filename for PPM format file
