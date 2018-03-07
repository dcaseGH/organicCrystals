''' Read a crystal optimiser z matrix file, and expand it for same molecule
    Useful when need to manually specify z matrix 
    Run as python [this script] [template res/cif file with high z'] [zmatrix file with one molecule in] > outputFile '''

import sys, re
#from relabellingString import labelIncludingNumberMols
from ccdc.io import CrystalReader

def labelIncludingNumberMols(atom1Label, stoichiometry, molIndex):
    element1, number1 = re.split('(\d+)', atom1Label)[:2]
    return element1 + str(int(number1) + int(stoichiometry[element1]) * molIndex)


blurb = '''This file contains the Z-matrix definition for each of the molecules
in the asymmetric unit of the crystal to be minimised by CrystalOptimizer.

It is absolutely necessary to contain an 'introductory' line for the Z-matrix
of each molecule followed imidiately (no empty line) by the Z-matrix
definition. The introductory line should be 'Z-matrix for molecule'
followed by the enumeration of teh molecule. The Z-matrix definition
should be done using the atom labels in exactly the same way as they
appear in the .res file.'''




crystal       = CrystalReader(sys.argv[1])[0]
stoichiometry = dict([re.split('(\d+)', x)[:2] for x in crystal.molecule.components[0].formula.split()])
zMatrixString = open(sys.argv[2], 'r').read()
print blurb#zMatrixString
for i in xrange(len(crystal.molecule.components)):
    print "\nZ-matrix for molecule %s"%(i + 1)
    for l in zMatrixString.split('\n')[1:]:
#        print l
        print "    ".join([labelIncludingNumberMols(x, stoichiometry, i) for x in l.split()])
