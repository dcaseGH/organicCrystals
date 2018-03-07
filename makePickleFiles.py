import re
import numpy as np
import pandas as pd
import pickle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pickle
from scipy.stats import gaussian_kde

from ase import Atoms
from ase.io import read as ASERead
from ase.spacegroup import crystal as ASECrystal
from ase.spacegroup import Spacegroup as ASESpacegroup

from ioAndInterfaces import ccdcCrystalToASE
from analyseClusters import generateNeighbouringAtoms, calculateClusterInformation, SetASEStructures
from analyseClusters import clusterInfoTo3DPoints, separateClusterInfo
#from ccdc.io import CrystalReader

dihedralIndices = {'FFA': [6, 17, 7, 8],
                   'MFA': [6, 15, 7, 8],
                   'TFA': [6, 14, 7, 8]}
groupLabels = {'FFA': {'ringA': [1,2,3,4,5,6,21,22,23,24],
                       'ringB': [7,8,9,10,11,12,26,27,28],
                       'ortho': [29],
                       'meta' : [13,14,15,16],
                       'NH'   : [17,25],
                       'COOH' : [0, 18, 19, 20]
                      },
               'TFA': {'ringA': [1,2,3,4,5,6,19,20,21,22],
                       'ringB': [7,8,9,10,11,12,24,25,26],
                       'ortho': [13,27,28,29],
                       'meta' : [17],
                       'NH'   : [14,23],
                       'COOH' : [0,15,16,18]},
               'MFA': {'ringA': [1,2,3,4,5,6,19,20,21,22],
                       'ringB': [7,8,9,10,11,12,24,25,26],
                       'ortho': [14,30,31,32],
                       'meta' : [13,27,28,29],
                       'NH'   : [15,23],
                       'COOH' : [0,16,17,18]}
              }
assert(sorted([x for k in groupLabels['FFA'].keys() for x in groupLabels['FFA'][k]]) == range(30))
assert(sorted([x for k in groupLabels['TFA'].keys() for x in groupLabels['TFA'][k]]) == range(30))
assert(sorted([x for k in groupLabels['MFA'].keys() for x in groupLabels['MFA'][k]]) == range(33))

mfaStructures =  SetASEStructures(ASERead('/Users/cmdc2/crystalStructureCoding/structureAnalysis/fenamatesClusterAnalysis/MFA_CSPPCM_137.cif', 
                                          index=':'),
                                  [1, 0 ,2], dihedralIndices['MFA'], groupLabels = groupLabels['MFA'])

print len(mfaStructures.structures)
tfaStructures =  SetASEStructures(ASERead('/Users/cmdc2/crystalStructureCoding/structureAnalysis/fenamatesClusterAnalysis/TFA_CSPPCM_137.cif', 
                                          index=':'),
                                  [1, 0 ,2], dihedralIndices['TFA'], groupLabels = groupLabels['TFA'])

print len(tfaStructures.structures)
ffaStructures =  SetASEStructures(ASERead('/Users/cmdc2/crystalStructureCoding/structureAnalysis/fenamatesClusterAnalysis/FFA_CSPPCM_137.cif', 
                                          index=':'),
                                  [1, 0 ,2], dihedralIndices['FFA'], groupLabels = groupLabels['FFA'])

print len(ffaStructures.structures)

pickle.dump(mfaStructures, open('mfaTemp.p', 'w'))
pickle.dump(tfaStructures, open('tfaTemp.p', 'w'))
pickle.dump(ffaStructures, open('ffaTemp.p', 'w'))
