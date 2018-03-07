usefulMessage = ''' Script to run as python [this script] [res file to be relabelled] [template res file] [output res name] - also add "y" as last arg to reorder and relabel things to standards '''

import sys
from ccdc.io import CrystalReader
from ccdc.descriptors import MolecularDescriptors
import re

def cpossPrecidence(label):
    ''' Turns atom label into a float that can be sorted
        Strictly should just make Carbons first, everything else, H last 
        Note potential bugs with tricky labels e.g. C1abc vs C1cab - last bit is untested '''

    #order is a bit arbitrary - chosen to try to fit current molecules (loosely)
    cpossOrder = {'C': 0.,
                  'F': 1000.,
                  'N': 2000.,
                  'O': 3000.,
                  'Cl': 4000.,
                  'H': 10000000.
                 }
    
    informationList = filter(None, re.split(r'(\d+)', label))
#    assert(len(informationList) < 3)
    if len(informationList) == 1:
        return cpossOrder[informationList[0]]
    elif len(informationList) == 2:
        return cpossOrder[informationList[0]] + float(informationList[1])
    #untested, not unique??, caveat
    elif len(informationList) == 3:
        return cpossOrder[informationList[0]] + float(informationList[1]) +\
            float(sum([ord(x)*100**(ix) for x in enumerate(informationList[2])])) *1.e-12

    raise Exception("label is dodgy - %s"%label)

def reorderLowestNumbering(listAtomData):
    ''' ensures numbers are ascending from 1 for each atom (maintains order) '''
    usedAtoms, modifiedAtomData = [], []
    for at in listAtomData:
        currentAtom = filter(None, re.split(r'(\d+)', at.split()[0].lower()))[0].upper()
#        if currentAtom in usedAtoms:
        modifiedAtomData.append(currentAtom + str(len([x for x in usedAtoms if x == currentAtom]) + 1) + \
                                "  " + "  ".join(at.split()[1:]))
        usedAtoms.append(currentAtom)
    
    return modifiedAtomData

def reorderedRes(oldRes,
                 _replacementDict,
                 enforceStandardListOrder = False,
                 enforceLowestNumbering = False):

    ''' Reads the oldRes, but changes labels according to _replacementDict 
        Returns a string
        enforeStandardListOrder means order is C1, C2.... N1, O1, O2, ... H1, H2 etc 
        enforceLowestNumbering drops numbers to be continuous from 1 (so if res originally
        started from 10, or 4 was missing -> 1,2,3,4,5 ... '''

    with open(oldRes, 'r') as inFile:
        oldData     = inFile.readlines()
        index       = [i for i, s in enumerate(oldData) if 'FVAR' in s][0]
        header      = "".join(oldData[:index+1])
        oldAtomData = "".join(oldData[index+1:])
        
        newAtomData = []
        for l in oldAtomData.split("\n"):
            
            if "END" in l:
                break
            
            tempL = l.split()
            tempL[0] = _replacementDict[tempL[0]]
        
            newAtomData.append(" ".join(tempL))

        if enforceStandardListOrder:
            newAtomData = sorted(newAtomData, key = lambda x: cpossPrecidence(x.split()[0]))
            
        if enforceLowestNumbering:
#            print newAtomData, enforceLowestNumbering, enforceLowestNumbering(newAtomData)
            newAtomData = reorderLowestNumbering(newAtomData)
            
        return header + "\n".join(newAtomData) + "\nEND\n"

def replacementDict(crystal1, crystal2):
    ''' Doesnt seem to consider 3d structure?? so may mix up enantiomers, or quasi-enantiomers
    ie those where the labelling of atoms makes them independent
    Presumably, as used as a dictionary, mol1 <-> mol2 makes little diff??? '''

    import re
    #only works with crystal1 being z'=1 and crystal2 being same molecule but many times
    assert(len(crystal1.molecule.components) == 1)
    assert(all([len(x.atoms) == len(crystal1.molecule.atoms) for x in crystal2.molecule.components]))
    
    def labelIncludingNumberMols(atom1Label, stoichiometry, molIndex):

        if len(re.split('(\d+)', atom1Label)) > 1:
            element1, number1 = re.split('(\d+)', atom1Label)[:2]

        # It may be possible that something doesn't have a number- give it '1' - watch here if other bugs arise
        else:
            element1, number1 = re.split('(\d+)', atom1Label)[0], u'1'

        return element1 + str(int(number1) + int(stoichiometry[element1]) * molIndex)

    atomsPerMol = dict([re.split('(\d+)', x)[:2] for x in crystal1.molecule.formula.split()])
     #    return dict([[labelIncludingNumberMols(pair[0].label, atomsPerMol, i2), pair[1].label]
    return dict([[pair[1].label, labelIncludingNumberMols(pair[0].label, atomsPerMol, i2)]
                 for mol1 in crystal1.molecule.components
                 for i2, mol2 in enumerate(crystal2.molecule.components)
                 for pair in MolecularDescriptors.MaximumCommonSubstructure().search(mol1, mol2)[0]])

def main():
    # names of files (all strings)
    try:
        try:
            oldResFile, templateResFile, newFilename = sys.argv[1:]
        except:
            oldResFile, templateResFile, newFilename, enforcement = sys.argv[1:]
    except:
        raise Exception(usefulMessage)

    # use the above subroutines to achieve the relabelling
    # uses ccdc molecule instances
    with open(newFilename, 'w') as outf:
        if len(sys.argv) == 4:
            outf.write(reorderedRes(oldResFile, replacementDict(CrystalReader(templateResFile)[0],
                                                                CrystalReader(oldResFile)[0]))
            )
        elif len(sys.argv) == 5 and sys.argv[4] == 'y':
            outf.write(reorderedRes(oldResFile, replacementDict(CrystalReader(templateResFile)[0],
                                                                CrystalReader(oldResFile)[0]),
                                    enforceStandardListOrder = True,
                                    enforceLowestNumbering   = True)
                       )
            
            
if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    main()

