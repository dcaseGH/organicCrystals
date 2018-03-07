# turn xyz into orient variables ASSUMING TRANSLATION ONLY
# also assume nanocrystal is written to xyz first, and surface is frozen

import numpy as np
from multipoleFile import BOHRTOANGSTROM
import sys

#assumed names (can give an xyz trajectory)
inputFile   = 'surfaceCalc.in'
if len(sys.argv) == 2:
    xyzFile     = sys.argv[1]
else:
    xyzFile     = 'output.xyz'
nanoFile    = 'nanoCryst.punch'
surfaceFile = 'surface.cell'


from ase.io import read as ASERead

firstAtomTrajectory = [image[0].position for image in ASERead(xyzFile, index=':')]

#quick hack is to read variables from input
with open(inputFile, 'r') as inFile:
    variables, startReading = [], False
    for l in inFile.readlines():
        if startReading and 'END' in l:
            break
        if l == 'VARIABLES\n':# in l:
            startReading = True
            continue
        if startReading:
            variables.append([l.split()[0], float(l.split()[1]), l.split()[2]])


assert(all([x[2]=='Bohr' for x in variables[:3]]))
initialVariables = np.array([x[1] for x in variables[:3]])
print 'Final translational variables (Bohr)', initialVariables - BOHRTOANGSTROM**-1 * (firstAtomTrajectory[0] - firstAtomTrajectory[-1])
