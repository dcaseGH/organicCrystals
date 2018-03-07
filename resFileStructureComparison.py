''' Run this script as python resFileStructureComparion.py [res file 1] [res file 2]
    If needed export PYTHONPATH='/home/dcase/davesScripts':$PYTHONPATH 
              export CSDHOME=/etc/CCDC/CSD_2017/ 
              export LD_LIBRARY_PATH=/usr/local/lib/python2.7/dist-packages/ccdc/_lib:$LD_LIBRARY_PATH '''

from ccdc.io import CrystalReader
from ccdc.crystal import PackingSimilarity

import sys

def returnNmatched_molecules(comparisonObject, crystal1, crystal2, returnRMSD=False):
    ''' All inputs are CCDC objects - this really just catches fails 
        Copied from davesScripts to save importing '''

    compareResult = comparisonObject.compare(crystal1, crystal2)

    if returnRMSD:
        try:
            return compareResult.nmatched_molecules, compareResult.rmsd
        except:
            return None, None
    else:
        try:
            return compareResult.nmatched_molecules
        except:
            return None
            
# Set up comparison object
ps = PackingSimilarity()
allowMolecularDifferences = True
clusterSize = 15
ps.settings.allow_molecular_differences = allowMolecularDifferences
ps.settings.packing_shell_size          = clusterSize

# define crystals
try:
    crystal1 =  CrystalReader(sys.argv[1], format='res')[0]
except:
    crystal1 =  CrystalReader(sys.argv[1], format='cif')[0]
try:
    crystal2 =  CrystalReader(sys.argv[2], format='res')[0]
except:
    crystal2 =  CrystalReader(sys.argv[2], format='cif')[0]

matchingData = returnNmatched_molecules(ps, crystal1, crystal2, returnRMSD=True)

#print "Matching {1} molecules (out of {2}) with RMSD = {3} Angstroms".format(clusterSize,clusterSize,clusterSize)
print "Allowing molecular differences? ", allowMolecularDifferences
print "Matching %s molecules (out of %s) with RMSD = %s Angstroms"%(matchingData[0],
                                                                    clusterSize,
                                                                    matchingData[1]) 
