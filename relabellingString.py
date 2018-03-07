usefulMessage = ''' Script to run as python [this script] [string to be relabelled] [old res file] [template res file] '''
import sys
import re
from ccdc.io import CrystalReader
from ccdc.descriptors import MolecularDescriptors
from relabellingRes import replacementDict

def relabelString(stringIn, replacementDict):
    ''' String is of the type entered into crystal optimizer '''
    #replacementDict is one to many - loop through elements
    print " ".join(["_".join(replacementDict[x].split("_")) for x in stringIn.split()])
    # names of files (all strings)

def labelIncludingNumberMols(atom1Label, stoichiometry, molIndex):

    element1, number1 = re.split('(\d+)', atom1Label)[:2]
    return element1 + str(int(number1) + int(stoichiometry[element1]) * molIndex)

def simpleRelabelling(stringIn, crystal):# replacementDict):
    ''' Simple relabelling assumes that you can just add nAtomType to each label
        i.e. if 2 O atoms you can assume O1 will be O3, O5, O7 in other molecules etc '''
#    atomsPerMol = dict([re.split('(\d+)', x)[:2] for x in crystal.molecule.formula.split()])
    atomsPerMol = dict([re.split('(\d+)', x)[:2] for x in crystal.molecule.components[0].formula.split()])
#    print atomsPerMol
#    print " ".join(["_".join(x.split("_")) for x in stringIn.split()])

#    def labelIncludingNumberMols(atom1Label, stoichiometry, molIndex):
#        element1, number1 = re.split('(\d+)', atom1Label)[:2]
#        return element1 + str(int(number1) + int(stoichiometry[element1]) * molIndex)
            
    for i2 in xrange(len(crystal.molecule.components)):
        print " +".join(["_".join([labelIncludingNumberMols(y, atomsPerMol, i2) for y in x.split("_")]) for x in stringIn.split()])
        
    return
            
        
    return dict([[pair[1].label, labelIncludingNumberMols(pair[0].label, atomsPerMol, i2)]
                 for mol1 in crystal1.molecule.components
                 for i2, mol2 in enumerate(crystal2.molecule.components)
                 for pair in MolecularDescriptors.MaximumCommonSubstructure().search(mol1, mol2)[0]])

def main():
    
    try:
#        oldString, oldResFile, templateResFile = sys.argv[1:]
        oldString, resFile = sys.argv[1:]
    except:
        raise Exception(usefulMessage)
#    print simpleRelabelling("+O1_C1_C2_C7 +H1_O1_C1_C2 +C1_C2_C7_N1 +C7_N1_C8_C9 +H6_N1_C8_C9 +F1_C14_C12_C11 +H6_N1_C8 +H1_O1_C1",
    print simpleRelabelling("O1_C1_C2_C7 H1_O1_C1_C2 C1_C2_C7_N1 C7_N1_C8_C9 H6_N1_C8_C9 F1_C14_C12_C11 H6_N1_C8 H1_O1_C1",
                            CrystalReader(resFile)[0])
if __name__ == '__main__':
    main()

    exit()
    print relabelString("+O1_C1_C2_C7 +H1_O1_C1_C2 +C1_C2_C7_N1 +C7_N1_C8_C9 +H6_N1_C8_C9 +F1_C14_C12_C11 +H6_N1_C8 +H1_O1_C1",
                        replacementDict(CrystalReader(templateResFile)[0],
                                        CrystalReader(oldResFile)[0]))
exit()


exit()

# use the above subroutines to achieve the relabelling
# uses ccdc molecule instances
with open(newFilename, 'w') as outf:
    outf.write(reorderedRes(oldResFile, replacementDict(CrystalReader(templateResFile)[0],
                                                        CrystalReader(oldResFile)[0]))
               )
