from ioAndInterfaces import ccdcCrystalToASE
from ccdc.io import EntryReader
from ase.io import write as ASEWrite

csdEntryReader = EntryReader('CSD')

aseCrystal = ccdcCrystalToASE(csdEntryReader.crystal('ABEBUF'))
ASEWrite('temp.xyz', aseCrystal)
