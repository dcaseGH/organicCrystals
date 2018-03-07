import numpy as np

class OrientInput():
    ''' General input for different Orient calculations '''
    def __init__(self):
        pass

    def allocateMemoryString(self):
        ''' Automate this when understand - dummy numbers for now '''
        return 'ALLOCATE\n'+\
               ' SITES %s\n'%(500)+\
               ' VARIABLES %s\n'%(5000)+\
               'END\n\n'
        
    def variablesString(self, variables):
        ''' Variables is a list, each element has two parts (variable name and value) '''
        if variables is None:
            return ''
        else:
            return 'VARIABLES\n' +\
                   '\n'.join([" %s  %s  %s"%(v[0], v[1], v[2]) for v in variables]) +\
                   '\nEND\n\n'

    def typesString(self, elements, addCustom = [['Hn', 'Z', 1]]):
        ''' Automate this '''
        from ase.data import atomic_numbers
        #custom H always Hn??
        dataList = [[e, 'Z', atomic_numbers[e]] for e in elements]
        dataList += addCustom
        #order?
        dataList = sorted(dataList, key=lambda x:x[2])
        print dataList
        return 'TYPES\n' +\
               '\n'.join([" %s  %s  %s"%(x[0], x[1], x[2]) for x in dataList]) +\
               '\nEND\n\n'

    def moleculeString(self, moleculeName, variableNames, moleculePunchFile):
        ''' Point to a punch file (could write it out in full) and explain orientation variables '''
        outString = "MOLECULE %s AT "%moleculeName + \
                    " ".join(["%s"%x for x in variableNames[:3] + ['ROTATED'] + variableNames[3:]]) + '\n\n'
        outString += "#include %s\n\n"%moleculePunchFile
        outString += 'END\n\n'
        return outString

    def nanocrystalString(self, cell, size):
        ''' This is just tiling of molecule in nanocrystal
            cell (np.array) is lower triangular - size is a x b x 1 (see assertions) '''

        assert(size[2] == 1)

        outString = 'NANOCRYSTAL\n'
        outString += self.cellString(cell)
        outString += '  SIZE A %s B %s C %s\nEND\n\n'%(size[0], size[1], size[2])
        return outString

    def cellString(self, cell):
        ''' Write a (lower triangular) cell '''
        assert(np.allclose(cell,
                           np.tril(cell)))
        outString  = '  CELL\n'
        # copy Rui's example way of writing cell
        outString += '    a   %s   %s\n'%(cell[0,0], cell[0,1])
        outString += '    b   %s   %s\n'%(cell[1,0], cell[1,1])
        outString += '    c   %s   %s   %s\n'%(cell[2,0], cell[2,1], cell[2,2])
        outString += '  END\n'
        return outString
        
    def includeFile(self, filename, message=None):
        if message:
            outString = message
        else:
            outString = ''
        return outString + '#include %s\n\n'%filename

    def sitesString(self, multipoleFileObjectList):
        outString = 'SITES\n'
        #must have given the object information to describe rotations
        for m in multipoleFileObjectList:
            assert (m.rotationInfo is not None)
            #set special atom types (acidic H) prior to this
            if any(['type' not in x.keys() or x['type'] is None for x in m.atomList]):
                m.setAtomTypes()
            outString += m.stringFormat(header = '', printTypes=True) + '\n'# header could be \n or END or something for z>1???
        outString += 'END\n\n'
        return outString

    def cellAndSitesString(self, aseCell, multipoleFileObjectList, header = '', direct = 60.0, reciprocal = 6):
        ''' give a cell (matrix) and a multipoleFileObject and concatenate them (basically) 
            direct and reciprocal are presumably Ewald instructions '''
        outString = header
        outString += 'LATTICE\n  RANGE DIRECT %s RECIPROCAL %s\n'%(direct, reciprocal)
        outString += self.cellString(aseCell)
        outString += self.sitesString(multipoleFileObjectList)
        outString += 'END'
        return outString

    
#def exampleNanocrystalInstructions(maxX, maxY, freezeRotation = False, surfaceXMin=-3, surfaceXMax=3, surfaceYMin=-3, surfaceYMax=3):
def exampleNanocrystalInstructions(maxX, maxY, surfaceMinA, surfaceMaxA, surfaceMinB, surfaceMaxB, freezeRotations = False):
    ''' maxX and maxY are the size of the nanocrystal which is written here so as to make trajectories '''
    outString = '''UNITS kj/mol Bohr

TIME

ENERGY

OPTIMIZE
  TITLE "OZPN 100 surface with a nanocrystal: Minimum"
  EVF VERBOSE
!  BFGS
  STEP 0.20 bohr
  SHOW steps
'''
    if freezeRotations:
        outString += '''  FIX
    dummyMoleculeName ORIENTATION
  END
'''
    outString += '''!  SHOW ENERGY VARIABLES
  ITERATIONS 240
!  CONVERGENCE 0.03
  PLOT XYZ FILE output.xyz EVERY 1 LATTICE %s %s %s %s NANOCRYSTAL 0 %s 0 %s
END

TIME

SHOW VARIABLES

FINISH'''%(surfaceMinA, surfaceMaxA, surfaceMinB, surfaceMaxB, maxX, maxY)
    return outString

# Points:
# variables should be precisely named
# why specify molecule, surface, cell etc - learn specific details
# cell settings
    
class NanocrystalOnSurfaceInput(OrientInput):
    ''' Optimise position of this nanocrystal '''
    def __init__(self):
        OrientInput.__init__(self)

# should have some way of writing the xyz to check that the surface is ok before running ??
#    def xyzString(self, 
#                  variable           = None,
#                  nano

    def inputString(self,
                    title            = 'This is a nanocrystal on surface calcn',
                    variables        = None,
                    nanocrystal      = None,
                    surfaceCrystal   = None,
                    nanocrystalSize  = None,
                    cellFile         = None,
                    potentialFile    = None,
                    moleculeFilename = None,
                    atomTypes        = None,
                    freezeRotations  = False
    ):

        ''' Possible issues - nanocrystal is ase atoms, surfaceCrystal is Crystal '''

        outString = 'NOTE ' + title + '\n\n'
        outString += self.allocateMemoryString()
        outString += self.variablesString(variables)
        outString += 'UNITS kj/mol Angstrom\n'
        outString += self.typesString(atomTypes)

        # The molecule is the nanocrystal - the nanocrystal details explain the tiling
        # (so put multiple molecules in molecule (if this makes sense))
        outString += self.moleculeString('dummyMoleculeName', [x[0] for x in variables], moleculeFilename)
        outString += 'NOTE give some details of nanocrystal and face\n\n'
        outString += self.nanocrystalString(nanocrystal.cell, nanocrystalSize)
        outString += self.includeFile(cellFile, message = 'NOTE include the lattice setup file here\n\n')
        outString += self.includeFile(potentialFile, message = 'NOTE include the pair potentials file here\n\n')

        # Delete example instructiosn when have more of an idea of what to do
        outString += exampleNanocrystalInstructions(*[nanocrystalSize[0] - 1,
                                                      nanocrystalSize[1] - 1]\
                                                    + self.boundsShowSurface(nanocrystalSize, 
                                                                             nanocrystal.cell, 
                                                                             surfaceCrystal.aseCell, 
                                                                             [x[1] for x in variables[3:6]]),
                                                     freezeRotations = freezeRotations)
        return outString

    def boundsShowSurface(self, nanocrystalSize, nanoCell, surfaceCell, angleVariables):
        ''' estimate bounds for surface during plotting (wont effect calculation) 
            returns surfaceMinA, surfaceMaxA, surfaceMinB, surfaceMaxB
            (status - run this and visually inspected things, but not extended to general angles, thoroughly tested) '''

        cornersNanoCrystal = np.dot(np.array([[nanocrystalSize[0], 0., 0.],
                                              [0., nanocrystalSize[1], 0.],
                                              [nanocrystalSize[0], nanocrystalSize[1], 0.],
                                              [0., 0., 0.]]),
                                    nanoCell)

        #zyz convention -- (after comparing to Stone's book) -- 
#        assert(all([x==0. for x in angleVariables[1:]]))
        from quaternions import zyz_quaternion, quaternion_rotatn
        cornersNanoCrystal = np.array([quaternion_rotatn(x, zyz_quaternion(*map(lambda x: x*np.pi/180., angleVariables))) for x in cornersNanoCrystal])

        #find bounds in surface cell
        boundsSurface = np.dot(cornersNanoCrystal,
                               np.linalg.inv(surfaceCell))
        boundsSurface = np.vstack([np.min(boundsSurface, axis=0),
                                   np.max(boundsSurface, axis=0)])

        return map(int, [np.floor(boundsSurface[0,0]),
                         np.ceil(boundsSurface[1,0]),
                         np.floor(boundsSurface[0,1]),
                         np.ceil(boundsSurface[1,1])])
      
