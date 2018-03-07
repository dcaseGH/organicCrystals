from ase.build import niggli_reduce
from ase import Atoms
from ase.io import write as aseWrite
from ioAndInterfaces import ccdcCrystalToASE
from ccdc.io import CrystalReader

#example input and outputs to test this
inputRes, outputRes = 'testingScripts/new.res', 'temp070917.res'

myCell = CrystalReader(inputRes)[0]
myASECell = ccdcCrystalToASE(myCell)
myASECell.set_pbc(True)
niggli_reduce(myASECell)

aseWrite(outputRes, myASECell)

print 'issue with this is that crystal optimiser wants labels etc, and probably fussy with res file format'
#could attach labels or something if can be bothered

print 'another issue is that the cell may not have the same convention, so would have to reset rather
       than just changing angles etc'
