''' Subroutines to print a list of closest molecules and information about rotations and
    close contacts etc '''

import numpy as np

# cell manipulation
def vectorLengths(matrix):
    return np.array([np.linalg.norm(v) for v in matrix])
        
def internalAngles(matrix):
    ''' '''
    _nRows = matrix.shape[0]
#    print _nRows, matrix, 
    return np.arccos(np.array([np.dot(matrix[(i+1)%_nRows], matrix[(i+2)%_nRows]) /\
                            (np.linalg.norm(matrix[(i+1)%_nRows]) *\
                             np.linalg.norm(matrix[(i+2)%_nRows])) for i in xrange(_nRows)]))

def latticeParameters(inCell):
    return list(vectorLengths(inCell)) + list((180./np.pi) * internalAngles(inCell))


def atomsToPDB(fileobj, atoms):
    ''' ASE idiotically truncates to 3dp when writing .pdb 
        This is a simplification of ase write (note - code taken from ase/io/proteindatabank.py
        Keep as atoms rather than Crystal to fit with aseToRDKit stuff '''
    from ase.geometry import cell_to_cellpar
    # THIS CANNOT BE READ BY RDKIT SO IS NOT MUCH USE !!!
    if isinstance(fileobj, basestring):
        fileobj = paropen(fileobj, 'w')

    cellpar = cell_to_cellpar(atoms.get_cell())
    fileobj.write('CRYST1 ' + ' '.join(['%s'%x for x in cellpar]) + ' P 1\n')
    fileobj.write('MODEL 1\n')

    symbols = atoms.get_chemical_symbols()
    natoms  = len(symbols)
    format = ('ATOM  %s %s MOL     1    %s %s %s 1.00 0.00'
              '          %2s  \n')
    p = atoms.get_positions()
    for a in range(natoms):
        x, y, z = p[a]
        fileobj.write(format % (a, symbols[a],
                                x, y, z, symbols[a].upper()))
    fileobj.write('ENDMDL\n')

    return fileobj

#move to ioAndInterfaces
def aseToRDKit(atoms, truncatedPositions = True):
    ''' Do via temp file with .pdb - note accuracy lost and formatting not good with atomsToPDB  '''

    from ase.io import write as ASEWrite
    from tempfile import SpooledTemporaryFile as stfile
    from rdkit import Chem

    #write a dummy file and read it to interface ASE and RDKit
    with stfile() as tempFile:
        if truncatedPositions:
            ASEWrite(tempFile, atoms, format='proteindatabank')
#        else:
#            tempFile = atomsToPDB(tempFile, atoms)
        tempFile.seek(0)
        return Chem.MolFromPDBBlock(tempFile.read(), sanitize=False, removeHs=False)

def rdKitToASE(rdKitMol):
    ''' Assume just one molecule (or do whole thing) '''

    from ase.io import read  as ASERead
    from tempfile import SpooledTemporaryFile as stfile
    from rdkit import Chem

    with stfile() as tempFile:
        tempFile.write(Chem.MolToPDBBlock(rdKitMol))
        tempFile.seek(0)
        return ASERead(tempFile, format='proteindatabank')

def manualSplitASEAtomsToMols(inAtoms, listAtomNumbers):
    ''' listAtomNumbers gives number atoms in each piece of ordered inAtoms list '''

    from ase import Atoms
    #list format is [[0:natoms1][natoms1+1:natoms1+natoms2]...]
    assert(all([x[0] < x[1] for x in listAtomNumbers]))
    assert(all([listAtomNumbers[i][1] == listAtomNumbers[i+1][0] for i in xrange(len(listAtomNumbers) - 1)]))

    outAtoms = []
    for x in listAtomNumbers:
        outAtoms.append(inAtoms[x[0]:x[1]])
        #keep the spacegroup and cell from inAtoms
        outAtoms[-1].info, outAtoms[-1].cell = inAtoms.info, inAtoms.cell
    return outAtoms

def splitASEAtomsToMols(inAtoms, useOriginalPositions=True):
    ''' Hack to turn one ASE atoms object into a list of molecules
        Definition of molecules comes from RDKit
        Always useOriginalPositions - ase writes low sig. fig. PDB files
        so save at start and copy back at end '''

    from ase.io import read  as ASERead
    from tempfile import SpooledTemporaryFile as stfile
    from rdkit import Chem

#    print 'inAtoms.get_positions()', inAtoms.get_positions()
    tempRDKitMols = aseToRDKit(inAtoms)

    #Do the opposite, generating molecules from GetMolFrags
    outAtoms = []
    for m in Chem.GetMolFrags(tempRDKitMols, asMols=True):
        with stfile() as tempFile:
            tempFile.write(Chem.MolToPDBBlock(m))
            tempFile.seek(0)
            outAtoms.append(ASERead(tempFile, format='proteindatabank'))
            #keep the spacegroup and cell from inAtoms
            outAtoms[-1].info, outAtoms[-1].cell = inAtoms.info, inAtoms.cell

    if useOriginalPositions:
        _tol = 1.e-2
#        _tol, _checkCounter = 1.e-2, 0
        for m in xrange(len(outAtoms)):
            for at in xrange(len(outAtoms[m])):
                for sp in inAtoms.get_positions():
#                    print m, at, sp, outAtoms[m][at].position
                    if np.linalg.norm(outAtoms[m][at].position - sp) < _tol:
#                        print outAtoms[m][at].position, sp, np.linalg.norm(outAtoms[m][at].position - sp), _tol, m, at
                        outAtoms[m][at].position = sp
#                        _checkCounter += 1
                        break
                else:
                    raise Exception('Could not match an atom position whilse using original position is splitASEAtomsToMols')
#    print '[x.get_positions() for x in outAtoms]', [x.get_positions() for x in outAtoms]
    return outAtoms

def getSpeciesListCIF(inputFile):
    ''' inputFile is the name of an input file                                                                                                    
        use ASE at this point to quickly parse file                                                                                               
        DOES NOTHING ABOUT FRACTIONAL OCCUPANCY '''

    from ase import Atoms
    from ase.spacegroup import Spacegroup as ASESpacegroup
    from ase.io.cif import parse_cif
    aseParser = parse_cif(inputFile)

    for name, c in aseParser:
        scaled_positions = np.array([c['_atom_site_fract_x'],
                                     c['_atom_site_fract_y'],
                                     c['_atom_site_fract_z']]).T

        atoms = Atoms(symbols = c['_atom_site_type_symbol'],
                      scaled_positions=scaled_positions,
                      pbc  = True,
                      cell = [c['_cell_length_a'], c['_cell_length_b'], c['_cell_length_c'],
                              c['_cell_angle_alpha'], c['_cell_angle_beta'], c['_cell_angle_gamma']],
                      info = {'spacegroup': ASESpacegroup(c['_symmetry_int_tables_number'])})

        yield atoms

class Molecule():
    ''' has aseAtoms and cluster environment information '''
    def __init__(self, aseAtoms):
        self.aseAtoms      = aseAtoms
        self.environment   = None
        self.listNeighbours = None
        self.axisSystem    = None

    #generalize this to include other possible types of interactions
    def indicesAcidicHydrogens(self, rXH = 1.3, searchedForBonds = ['O', 'N']):
        ''' return indices of self.aseAtoms which are H and close (<rXH) to O or N '''
        indices = []
        bondDistances = self.aseAtoms.get_all_distances(mic=False)
        for ia, a in enumerate(self.aseAtoms):
            if a.symbol == 'H':
                if self.aseAtoms[np.where(bondDistances[:, ia] < rXH)[0][0]].symbol in searchedForBonds:
                    indices.append(ia)
        return indices

    def rotateASEAtomsAroundCentroidQuat(self, quat):
        from quaternions import quaternion_rotatn
        oldPositions = self.aseAtoms.get_positions()
        _mean = np.mean(oldPositions, axis=0)
        oldPositions -= _mean
        newPositions = [quaternion_rotatn(x, quat) for x in oldPositions]

        for i in xrange(len(self.aseAtoms)):
            self.aseAtoms[i].position = newPositions[i] + _mean

    def rotateASEAtoms(self, _R):#, moveMolecules = False):
        ''' Make sure _R has determinant of 1. (mirror ok if wanted I suppose) '''
        for a in self.aseAtoms:
            a.position = np.dot(a.position, _R)

    def translateASEAtoms(self, v):
        ''' cartesian translation '''
        for a in self.aseAtoms:
            a.position += v

    # the stuff below relates to the pictures of the atoms and local environments
    def calculateAxisSystem(self, atom1, atom2, atom3):
        ''' atom1 (int) is (index for) pivot, to atom2 is x, 3 -> y . Make right handed frame '''

        from quaternions import make_quat_rot2, quaternion_rotatn

        internalX = (self.aseAtoms[atom2].position - self.aseAtoms[atom1].position) /\
                     np.linalg.norm(self.aseAtoms[atom1].position - self.aseAtoms[atom2].position)
        internalY = (self.aseAtoms[atom3].position - self.aseAtoms[atom1].position) /\
                     np.linalg.norm(self.aseAtoms[atom3].position - self.aseAtoms[atom1].position)
        internalZ = np.cross(internalX, internalY) / np.linalg.norm(np.cross(internalX, internalY))
        _angle = np.pi/2. - vecAngle(internalX,
                                     internalY)
        rotationQ = make_quat_rot2(_angle, internalZ)
        internalY = quaternion_rotatn(internalY,
                                      rotationQ)

        self.axisSystem = np.vstack([internalX, internalY, internalZ])
        return self.axisSystem

    def calculateEnvironment(self, dihedralAtoms, groupLabels):
        ''' Assuming neighbours and axisSystem set already (+of course aseAtoms) '''

        from quaternions import inverse_quat_rot2Safe 

        outDict = {'dihedral': [],
                   'r(com)'  : [],
                   'v(com)'  : [],
                   'r(contact)':[],
                   'groups': [],
                   'detOperation': [],
                   'thetaRot': [],
                   'axisRot': []}

        #need to have defined a RH axis system metric 1^3
        assert(abs(np.linalg.det(self.axisSystem) - 1.) < 1.e-6)
        outDict['originalDihedral'] = self.aseAtoms.get_dihedral(*dihedralAtoms)
        
        _com = self.aseAtoms.get_center_of_mass()
        for x in self.listNeighbours:
            dihedralAngle = x.get_dihedral(*dihedralAtoms)
            comDistance   = np.linalg.norm(x.get_center_of_mass() - _com)
            comVector     = np.dot(self.axisSystem,
                               x.get_center_of_mass() - _com)
            closestContactInfo = closestContact(self.aseAtoms,
                                                x)
            if groupLabels:
#                adjacentGroups = ([k for k in groupLabels.keys() if closestContactInfo[1][0][0]%33 in groupLabels[k]][0],\
#                                  [k for k in groupLabels.keys() if closestContactInfo[1][1][0]%33 in groupLabels[k]][0])
                adjacentGroups = ([k for k in groupLabels.keys() if closestContactInfo[1][0][0] in groupLabels[k]][0],\
                                  [k for k in groupLabels.keys() if closestContactInfo[1][1][0] in groupLabels[k]][0])

            rotationData = rotationAnalysis(self.aseAtoms, x)
            _theta, _axis = inverse_quat_rot2Safe(rotationData[2])
#            assert(rotationData[1] < 1.e-3) #otherwise z'=2 ???? #want z'>1 crystals

            outDict['detOperation'].append(rotationData[0])
            outDict['thetaRot'].append(_theta)
#            outDict['axisRot'].append(_axis)
            outDict['axisRot'].append(np.dot(self.axisSystem,
                                             _axis))
            outDict['dihedral'].append(dihedralAngle)
            outDict['r(com)'].append(comDistance)
            outDict['v(com)'].append(comVector)
            outDict['r(contact)'].append(closestContactInfo[0])
            outDict['groups'].append(adjacentGroups)

        #make the lists of numbers/arrays into np.array
        for k in ['detOperation', 'thetaRot', 'axisRot', 'dihedral', 'r(com)', 'v(com)', 'r(contact)']:
            outDict[k] = np.array(outDict[k])
        self.environment = outDict
        return self.environment

    def environmentPointList(self):
        ''' Returns a list of points with environment '''
        return [dict([(k, self.environment[k][i]) for k in self.environment.keys() if k != 'originalDihedral']) 
                                                      for i in xrange(len(self.environment['thetaRot']))]
    def writeFile(self, filename):     
        ''' Write molecule to P1 cell to check things -note that ase cif writer might wrap atoms '''
        from ase import Atoms
        from ase.spacegroup import Spacegroup
        dummyAtoms = Atoms(symbols = self.aseAtoms.get_chemical_symbols(),
                           positions = self.aseAtoms.get_positions(),
                           cell      = 50. * np.eye(3),
                           info = {'spacegroup': Spacegroup(1)})
        dummyAtoms.write(filename)


class Crystal():
    ''' Asymmetric unit contains aseCrystals 
        Keep a cell/spacegroup outside if needed'''

    #should really be able to init another way...
    def __init__(self, aseAtomsList):
        self.asymmetricMolecules = [Molecule(x) for x in aseAtomsList]
        self.aseCell = aseAtomsList[0].cell
        self.aseInfo = aseAtomsList[0].info
        # also have something to hold calculated info e.g. pcmData.energyPerMol etc


    @classmethod
    def fromRes(cls, resFilename, spacegroup, reduceToP1=True):
        ''' Return a P1 cell '''
        from ioAndInterfaces import resToASEAtoms

        if not reduceToP1:
            raise Exception('Not implemented this yet')

#        print 'list(getSpeciesListCIF(cifFilename))', list(getSpeciesListCIF(cifFilename)), list(getSpeciesListCIF(cifFilename))[0].get_positions()
        _crystal = cls([resToASEAtoms(resFilename, spacegroup)])
#        print 'inside fromCif ',len(_crystal.asymmetricMolecules)
        if reduceToP1:
            _crystal.asymmetricMolecules = _crystal.filledUnitCellMolecules()

        return _crystal

    @classmethod
    def fromCif(cls, cifFilename, reduceToP1=True):
        ''' Return a P1 cell '''

        if not reduceToP1:
            raise Exception('Not implemented this yet')

#        print 'list(getSpeciesListCIF(cifFilename))', list(getSpeciesListCIF(cifFilename)), list(getSpeciesListCIF(cifFilename))[0].get_positions()
        _crystal = cls(splitASEAtomsToMols(list(getSpeciesListCIF(cifFilename))[0]))
#        print 'inside fromCif ',len(_crystal.asymmetricMolecules)
        if reduceToP1:
            _crystal.asymmetricMolecules = _crystal.filledUnitCellMolecules()

        return _crystal

    @classmethod
    def fromZp1Cif(cls, cifFilename, fillCell=False, manuallyBreakMoleculeList=None, nAtomsPerMol=None):
        ''' If you have a cif with a z'=1 structure, you can use this
            Otherwise need to use rdKit or something to break it up or
            give manuallyBreakMoleculeList to specify this (so cif isnt really z'=1...) '''

        _crystal = cls([list(getSpeciesListCIF(cifFilename))[0]])
        
        #dont do fill cell and manuallyBreak as untested
        if manuallyBreakMoleculeList is not None:
            _crystal.splitAsymmetricMoleculesManually(manuallyBreakMoleculeList)
        elif nAtomsPerMol is not None:
            nMols = len(_crystal.asymmetricMolecules[0].aseAtoms) / nAtomsPerMol
            _crystal.splitAsymmetricMoleculesManually([[i * nAtomsPerMol, (i+1)*nAtomsPerMol] for i in xrange(nMols)])
        elif fillCell:
            _crystal.asymmetricMolecules = _crystal.filledUnitCellMolecules()

        return _crystal
        
    def uniqueElements(self):
        ''' Whats it made of '''
        return set([x.symbol for m in self.asymmetricMolecules for x in m.aseAtoms])

    def generateNeighbouringAtoms(self, nMolecules = 15):
        ''' Generate all neighbouring molecules within certain unit cells and return
            a list of nMolecules which are closest by centre of mass distance 
            Atoms in name refers to molecules (ie ase Atoms) '''

        # can optimize speed by choosing translations at a later date if needed

        from ase import Atoms
        listNeighbours = []
        for i1, asymMol in enumerate(self.asymmetricMolecules):
            moreAtoms = []
            for i2, asymNeighbour in enumerate(self.asymmetricMolecules):
                moreAtoms += [Atoms(symbols = asymNeighbour.aseAtoms.get_chemical_symbols(),
                              scaled_positions = np.array([np.dot(rot, x) + trans for x in asymNeighbour.aseAtoms.get_scaled_positions(wrap=False)]),
                              cell=asymMol.aseAtoms.cell)
                              for rot, trans in extendTranslations(asymMol.aseAtoms.info['spacegroup'].get_symop())]

            self.asymmetricMolecules[i1].listNeighbours = sorted(moreAtoms,
                      key=lambda x: np.linalg.norm(x.get_center_of_mass() - asymMol.aseAtoms.get_center_of_mass()))[1:nMolecules+1]

        return

    def calculateAxisSystems(self, axisAtoms):
        for m in self.asymmetricMolecules:
            m.calculateAxisSystem(*axisAtoms)

    def calculateClusterInformation(self, axisAtoms, dihedralAtoms, groupLabels = None, dihedralMask = None, recalcNeighbours = False):
        ''' Makes a dictionary describing local environment of each asym mol 
            For now, axisAtoms/dihedralAtoms etc same for each mol '''

        if not all([x.listNeighbours is not None for x in self.asymmetricMolecules]) or recalcNeighbours:
            self.generateNeighbouringAtoms()#nMol = 15??

        if not all([x.axisSystem is not None for x in self.asymmetricMolecules]) or recalcNeighbours:
            self.calculateAxisSystems(axisAtoms)#axisSystem = molecularAxisSystem(inAtoms, *axisAtoms)

        for m in self.asymmetricMolecules:
            m.calculateEnvironment(dihedralAtoms, groupLabels)

        return

    def listAllMolecularEnvironments(self):
        ''' List of every environment point of every asym mol in every structure '''

        return [y for s in self.structures
                  for x in s['crystal'].asymmetricMolecules
                  for y in x.environmentPointList()]

    def filledUnitCellMolecules(self, superCell = np.ones(3), cellLimits = None):
        ''' Returns list of molecules that would fill supercell/bounds
            If cellLimits not set, use supercell (superflouous really) - must give one or just default unit cell '''

        from ase import Atoms

        outAtoms = []
        #take care with cell limits - typically for molecule [0, 1)^3 matrix part can make negative, but vec (trans) can increase 
        if cellLimits is None:
            cellLimits = np.vstack([np.array([superCell[0], superCell[1], superCell[2]]),
                                    np.zeros(3)])

        # cellLimits must be max over min, 3,2 shape
        assert(cellLimits.shape == (2,3))
        assert(all([cellLimits[0, i] > cellLimits[1, i] for i in xrange(3)]))

        #cell bounds are rounded up or down ints
        cellBounds = np.vstack([1 + np.ceil(cellLimits[0,:]),
                               -1 + np.floor(cellLimits[1,:])]).astype(int)

        for i1, asymMol in enumerate(self.asymmetricMolecules):
            outAtoms += [Atoms(symbols = asymMol.aseAtoms.get_chemical_symbols(),
                               scaled_positions = np.array([np.dot(rot, x) + trans for x in asymMol.aseAtoms.get_scaled_positions(wrap=False)]),
                               info=asymMol.aseAtoms.info,
                               cell=asymMol.aseAtoms.cell)
                         for rot, trans in extendTranslations(asymMol.aseAtoms.info['spacegroup'].get_symop(),
                                                              xLimits= cellBounds[:, 0],
                                                              yLimits= cellBounds[:, 1],
                                                              zLimits= cellBounds[:, 2])]
#            outAtoms = [x for x in outAtoms if (x.get_center_of_mass(scaled=True) >= cellLimits[1,:]).all() and\
#                        (x.get_center_of_mass(scaled=True) < cellLimits[0,:]).all()]

        outAtoms = [x for x in outAtoms if (x.get_center_of_mass(scaled=True) >= cellLimits[1,:]).all() and\
                                           (x.get_center_of_mass(scaled=True) < cellLimits[0,:]).all()]

        #check that found all z * z' * no. cells molecules 
        # this may not be true if not integer cell limits
        if (np.prod(cellLimits[0, :] - cellLimits[1, :])).is_integer():
            assert(len(self.asymmetricMolecules)*\
                   len(self.asymmetricMolecules[0].aseAtoms.info['spacegroup'].get_symop())*\
                   np.prod(cellLimits[0, :] - cellLimits[1, :]) == len(outAtoms))
        else:
            print 'not integer', np.prod(cellLimits[0, :] - cellLimits[1, :]), np.prod(cellLimits[0, :] - cellLimits[1, :]).is_integer()
            pass # could check its within bounds if ever feel the need
        return [Molecule(x) for x in outAtoms]

    def writeASEAtoms(self, filename, pbc=True, wrapAtoms=False):
        ''' Combine aseAtoms objects and then write to filename 
            N.B. Atoms.write puts things into P1 for some reason '''
        from ase import Atoms

        if filename.lower()[-4:] == '.xyz':
            from ase.spacegroup import Spacegroup
            dummyAtoms = Atoms(symbols   = np.concatenate([x.aseAtoms.get_chemical_symbols() for x in self.asymmetricMolecules]),
                               positions = np.concatenate([x.aseAtoms.get_positions() for x in self.asymmetricMolecules]),
                               cell      = 50. * np.eye(3),
                               info      = {'spacegroup': Spacegroup(1)})
            dummyAtoms.write(filename)
            return True

        _writingAtoms = Atoms(symbols = [x for m in self.asymmetricMolecules for x in m.aseAtoms.get_chemical_symbols()],
                              scaled_positions = np.array([x for m in self.asymmetricMolecules for x in m.aseAtoms.get_scaled_positions(wrap=False)]),
                              cell    = self.aseCell,
                              pbc     = pbc,
                              info    = self.aseInfo)
#                              cell    = self.asymmetricMolecules[0].aseAtoms.cell,
#                              pbc     = True,
#                              info    = self.asymmetricMolecules[0].aseAtoms.info)
        if wrapAtoms == False:
            from ioAndInterfaces import aseWriteCifNoWrap
            aseWriteCifNoWrap(filename, [_writingAtoms])
        else:
            _writingAtoms.write(filename)

    def setASECell(self, cell, changeMolecules = True):
        ''' Safe way to set cell inc for all molecules '''
        self.aseCell = cell
        for m in self.asymmetricMolecules:
            m.aseAtoms.cell = cell

    def setASEInfo(self, info, changeMolecules = True):
        ''' Safe way to set info (i.e. spacegroup) inc for all molecules '''
        self.aseInfo = info
        for m in self.asymmetricMolecules:
            m.aseAtoms.info = info

    def rotateASEAtoms(self, _R, moveMolecules = False):
        ''' Rotate all atoms by _R (and move back into unit cell if required) '''

        for m in self.asymmetricMolecules:
            m.rotateASEAtoms(_R)

    def translateASEAtoms(self, vec):
        ''' Translate all atoms by vec '''
        for m in self.asymmetricMolecules:
            m.translateASEAtoms(vec)

    def putAllMoleculesInUnitCell(self):
        ''' Moves molecules whose centre of mass isnt [0,1)^3 '''
        for m in self.asymmetricMolecules:
            m.translateASEAtoms(np.dot(np.mod(m.aseAtoms.get_center_of_mass(scaled=True), 1.) -\
                                       m.aseAtoms.get_center_of_mass(scaled=True),
                                       self.aseCell))

    def splitAsymmetricMoleculesManually(self, listAtomNumbers):
        ''' From a fake z'=1 structure and a list which specifies blocks into each
            new molecule, break the asymmetric molecule up '''

        self.asymmetricMolecules = [Molecule(x) for x in  manualSplitASEAtomsToMols(self.asymmetricMolecules[0].aseAtoms,
                                                                                    listAtomNumbers)]


    #def alignCell()
    def clone(self, matrix):
        ''' Give a matrix which can rotate axes amongst each other, then make a new crystal in this setting '''

        from ase import Atoms

        #only works for P1 (assert info['spacegroup'] == 1 (or something)
        # |matrix| = 1
        assert(abs(np.linalg.det(matrix) - 1.) < 1.e-7)

        _newLatticeParameters = latticeParameters(np.dot(matrix,
                                                         self.aseCell))

        print self.asymmetricMolecules[0].aseAtoms.get_scaled_positions(wrap=False)[0], np.dot(matrix,
                                                                                               self.asymmetricMolecules[0].aseAtoms.get_scaled_positions(wrap=False)[0])
        print np.dot(matrix,
                     self.asymmetricMolecules[0].aseAtoms.get_scaled_positions(wrap=False).T).T[0]
        _newCrystal = Crystal([Atoms(symbols   = m.aseAtoms.get_chemical_symbols(),
                                     scaled_positions = np.dot(matrix,
                                                        m.aseAtoms.get_scaled_positions(wrap=False).T).T,
                                     cell      = _newLatticeParameters) for m in self.asymmetricMolecules])

        return _newCrystal

class Surface(Crystal):
    ''' Surface is a crystal with normal vector etc '''
    def __init__(self, aseAtomsList, millerIndex = None):
        Crystal.__init__(self, aseAtomsList)
        self.millerIndex = millerIndex
        if self.millerIndex is not None:
            self.createCut()
#        print type(self)

#class method to generate a new instance (with all the molecules in, and keep old thing if needed??)
    def createCut(self):
        ''' if you give a fractional description, it will keep pbc (does Orient use these??) in 2D and a gap
            alternatively can just give a finite chunk '''
        from ase.spacegroup import Spacegroup as ASESpacegroup
        #ase has a reduce thing to stop miller indices getting too high if wanted
        #define surfaceBravaisLattice (the one with v1, v2 in the plane and v3 || millerIndex)
        self.surfaceBravaisLattice = self.surfaceDefinedLattice()

#        self.rotateASEAtoms(np.dot(np.linalg.inv(self.aseCell), self.surfaceBravaisLattice))
        #should really get a reduced surface (but note that only really makes sense when all three vectors defined)
        # if change cell, need to rotate molecules??? or can just get away without??

        #make enough molecules to fill cell
        # cell limits should be supercell of slab proj into inv cell
        _cellLimits = np.vstack([np.dot(np.dot(np.ones(3),
                                               self.surfaceBravaisLattice),
                                        np.linalg.inv(self.aseCell)),
                                 np.zeros(3)])
        _cellLimits = np.array([[max(_cellLimits[0,i], _cellLimits[1,i]) for i in xrange(3)],
                                [min(_cellLimits[0,i], _cellLimits[1,i]) for i in xrange(3)]])
        #at some point extend cell in 3rd vector or bravais lattice somehow
#        self.asymmetricMolecules = self.filledUnitCellMolecules(cellLimits = np.array([[4., 2., 1.], [0. ,0. ,0.]]))#[-1.e-6, -1.e-6, -1.e-6]]))
        #print len(self.filledUnitCellMolecules()), 1./( np.linalg.det(self.aseCell) / np.linalg.det(self.surfaceBravaisLattice)), len(self.filledUnitCellMolecules(cellLimits = _cellLimits))
        self.asymmetricMolecules = self.filledUnitCellMolecules(cellLimits = _cellLimits)
#        self.rotateASEAtoms(np.linalg.inv(self.aseCell), self.surfaceBravaisLattice) 
        #print 'new com frac ', np.dot([x.aseAtoms.get_center_of_mass() for x in self.asymmetricMolecules], np.linalg.inv(self.surfaceBravaisLattice))
        self.setASECell(self.surfaceBravaisLattice)
        print 'new com frac ', np.dot([x.aseAtoms.get_center_of_mass() for x in self.asymmetricMolecules], np.linalg.inv(self.aseCell))
        self.putAllMoleculesInUnitCell()
        print 'new com frac ', np.dot([x.aseAtoms.get_center_of_mass() for x in self.asymmetricMolecules], np.linalg.inv(self.aseCell))
#        self.asymmetricMolecules = self.filledUnitCellMolecules()
#        print len(self.asymmetricMolecules)
        #break symmetry (make P1)
        self.setASEInfo({'spacegroup': ASESpacegroup(1)})
#        if need to match up to give 2D
        # if doing cartesian chunk...
        #make cell very big 
#        self.setASECell(np.ones(3) * 50.)
        
        pass

    #could be static method if give a cell too??
    def surfaceDefinedLattice(self, millerIndex=None):
        ''' Uses equations from http://dx.doi.org/10.1016/j.susc.2013.05.016 to define
            a lattice with v1 and v2 in surface, and v3 (roughly) orthogonal (see paper above) '''

        # Use one already defined if not explicitly set
        if millerIndex is None:
            millerIndex = self.millerIndex

        def lcm(a, b):
            from fractions import gcd
            """Return lowest common multiple. From https://stackoverflow.com/questions/147515/least-common-multiple-for-3-or-more-numbers"""
            #need integer - take care with // as python 2 vs 3 issues
            return a * b // gcd(a, b)

        #case with no 0s in millerIndex
        if np.prod(millerIndex) != 0:
            print millerIndex
            _M  = float(reduce(lcm, list(millerIndex)))
            _points = np.vstack([(_M / float(millerIndex[0])) * self.aseCell[0],
                                 (_M / float(millerIndex[1])) * self.aseCell[1],
                                 (_M / float(millerIndex[2])) * self.aseCell[2]])

        #case with 1 0 in millerIndex
        elif len([x for x in millerIndex if x==0]) == 1:
            _zeroIndex = np.where(millerIndex == 0)[0][0]
            _M = float(lcm(millerIndex[(_zeroIndex + 1) % 3],
                           millerIndex[(_zeroIndex + 2) % 3]))
            _points = np.zeros((3,3))
            _points[0] = (_M / float(millerIndex[(_zeroIndex + 1) % 3])) * self.aseCell[(_zeroIndex + 1) % 3]
            _points[1] = (_M / float(millerIndex[(_zeroIndex + 2) % 3])) * self.aseCell[(_zeroIndex + 2) % 3]
            _points[2] = _points[0] + self.aseCell[_zeroIndex]
#            _points[(_zeroIndex + 1) % 3] = (_M / float(millerIndex[(_zeroIndex + 1) % 3])) * self.aseCell[(_zeroIndex + 1) % 3]
#            _points[(_zeroIndex + 2) % 3] = (_M / float(millerIndex[(_zeroIndex + 2) % 3])) * self.aseCell[(_zeroIndex + 2) % 3]
#            _points[_zeroIndex] = _points[(_zeroIndex + 1) % 3] + self.aseCell[millerIndex[_zeroIndex]]
#            print millerIndex, _zeroIndex, _M, _points
        #case with 2 0s in millerIndex
        elif len([x for x in millerIndex if x!=0]) == 1:
            _nonZeroIndex = np.where(millerIndex != 0)[0][0]
            _points = np.zeros((3,3))
            #points [ nonZeroIndex ] = 000
            _points[1] = self.aseCell[(_nonZeroIndex + 1) % 3]
            _points[2] = self.aseCell[(_nonZeroIndex + 2) % 3]
#            _points[(_nonZeroIndex + 1) % 3] = self.aseCell[(_nonZeroIndex + 1) % 3]
#            _points[(_nonZeroIndex + 2) % 3] = self.aseCell[(_nonZeroIndex + 2) % 3]
            #print _nonZeroIndex, _points, self.aseCell 
        #the bestNormal vector is most orthogonal to the surface
        _searchLim = 6 #only look in subset of possible vectors
        vectorList = [np.dot(np.array([i,j,k]), self.aseCell) for i in xrange(-_searchLim, _searchLim+1)
                                                              for j in xrange(-_searchLim, _searchLim+1)
                                                              for k in xrange(-_searchLim, _searchLim+1) if not (i==0 and j==0 and k==0)]

        absoluteNormal = np.cross(_points[1] - _points[0], _points[2] - _points[0])
        absoluteNormal /= np.linalg.norm(absoluteNormal)
        originalSign = np.sign(np.linalg.det(self.aseCell))
        bestNormal, bestValue = np.inf * np.ones(3), -1. * originalSign * np.inf

        _delta = 1.e-7
        for v in vectorList:
            currentLength = np.linalg.norm(v)
            currentValue  = originalSign * np.dot(v, absoluteNormal)/currentLength
            # if almost same value, take shorter vector (i.e. [100] better than [200])
            if currentValue > bestValue + _delta or\
               abs(currentValue - bestValue) < _delta and currentLength < np.linalg.norm(bestNormal):
                bestNormal, bestValue = v, currentValue 
        #maintain sign of coordinate system

        return np.vstack([_points[1] - _points[0], _points[2] - _points[0], bestNormal])
                                             
def closestContact(mol1, mol2):
    ''' Input two molecules (ase atoms)
        return (closest contact, array(atom in 1, atom in 2 that its close to) '''

    #should have van der waals radii
    distanceMatrix = np.array([[np.linalg.norm(at1.position - at2.position) for at1 in mol1] for at2 in mol2])
    closestPoint = np.min(distanceMatrix)
    return closestPoint, np.where(distanceMatrix == closestPoint)

def pureTranslations(xLimits=[-1,1],
                     yLimits=[-1,1],
                     zLimits=[-1,1]):
    return [(np.eye(3), np.array([x,y,z])) for x in xrange(min(xLimits), max(xLimits)+1)
            for y in xrange(min(yLimits), max(yLimits)+1)
            for z in xrange(min(zLimits), max(zLimits)+1)]

def extendTranslations(originalOperations,
                       xLimits=[-3,3],
                       yLimits=[-3,3],
                       zLimits=[-3,3]):
    #move to Crystal
    return [(originalOperation[0], originalOperation[1] + np.array([x,y,z]))
            for x in xrange(min(xLimits), max(xLimits)+1)
            for y in xrange(min(yLimits), max(yLimits)+1)
            for z in xrange(min(zLimits), max(zLimits)+1)
            for originalOperation in originalOperations]

def vecAngle(v1, v2):
    ''' v1-v2 angle'''
    return np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))

def rotationAnalysis(m1, m2, dihedral=None):
    ''' Rotation of 2 onto 1 with dihedral change if required 
        m1 and m2 are ase atoms '''
    from listMathsSlow import overlay_points_RMSD
    from copy import deepcopy
    _m2 = deepcopy(m2)

    rmsdInfo = overlay_points_RMSD(np.array([x.position for x in m1])  - np.mean(np.array([x.position for x in m1]), axis=0),
                                   np.array([x.position for x in _m2]) - np.mean(np.array([x.position for x in _m2]), axis=0))
    #z' == 1 !!!!! Must be the case here - ie no differences in internal angles
    if rmsdInfo[0] < 1.e-3:
        return (1., rmsdInfo[0], rmsdInfo[1])

    rmsdInfo = overlay_points_RMSD(np.array([x.position for x in m1])  - np.mean(np.array([x.position for x in m1]), axis=0),
                                   -1. *( np.array([x.position for x in _m2]) - np.mean(np.array([x.position for x in _m2]), axis=0)))
    return (-1., rmsdInfo[0], rmsdInfo[1])


def getMarker(info):
    if info == 1.:
        return 'x'
    elif info == -1.:
        return 'o'
    
def getColour(axis, theta):
    return tuple(np.clip(map(lambda x: abs(theta * x / np.pi), axis), 0., 1.))
    
    
def clusterInfoTo3DPointsOld(ax, clusterInfo):

    #temp subroutine until incorporated into something below
    for pt in xrange(len(clusterInfo['v(com)'])):

        ax.scatter(clusterInfo['v(com)'][pt][0],
                   clusterInfo['v(com)'][pt][1],
                   clusterInfo['v(com)'][pt][2],
                   c= getColour(clusterInfo['axisRot'][pt],
                                clusterInfo['thetaRot'][pt]),
                   marker=getMarker(clusterInfo['detOperation'][pt])) 

    return ax

def clusterInfoTo3DPoints(ax, listPoints):

    #temp subroutine until incorporated into something below
    for pt in listPoints:

        ax.scatter(pt['v(com)'][0],
                   pt['v(com)'][1],
                   pt['v(com)'][2],
                   c=getColour(pt['axisRot'],
                               pt['thetaRot']),
                   marker=getMarker(pt['detOperation'])) 

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    return ax

def clusterInfoTo3DPointsDihedral(ax, listPoints):
    ''' colour by dihedral '''

    #position and colour defined
    _x, _c = np.array([pt['v(com)']   for pt in listPoints]),\
             np.array([pt['dihedral'] for pt in listPoints])
    _p = ax.scatter(_x[:,0], _x[:,1], _x[:,2], c=_c)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
#    colorbar()
#    ax.legend()
    return _p

def clusterInfoTo3DPointsRotAxis(ax, listPoints):
    ''' colour by rotAxis '''

    #position and colour defined
    _x, _c = np.array([pt['v(com)']   for pt in listPoints]),\
             np.array([getColour(pt['axisRot'],
                               pt['thetaRot']) for pt in listPoints])
    _p = ax.scatter(_x[:,0], _x[:,1], _x[:,2], c=_c)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
#    colorbar()
#    ax.legend()
    return _p

#def separateClusterInfo(clusterInfo):
#    ''' Separate into detOperation, theta = 0, and pi 1,0,0 , 0,1,0 , 0,0,1 rotation (8 in tot) '''
#    def getColour(axis, theta):
#        return tuple(np.clip(map(lambda x: abs(theta * x / np.pi), axis), 0., 1.))
#
#    #temp subroutine until incorporated into something below
#    for pt in xrange(len(clusterInfo['v(com)'])):
#        ax.scatter(clusterInfo['v(com)'][pt][0],
#                   clusterInfo['v(com)'][pt][1],
#                   clusterInfo['v(com)'][pt][2],
#                   c= getColour(clusterInfo['axisRot'][pt],
#                                clusterInfo['thetaRot'][pt]),
#                   marker=getMarker(clusterInfo['detOperation'][pt])) 

#    return ax

#class Cluster():
#    ''' Holds information and subroutines for calculating and plotting clusters '''
#    def __init__(self):

# rename - not aseStructures - thats input - 
class SetASEStructures():
    ''' Holds set of ASE structures - extract structural data '''

    def __init__(self, structureFilename, axisAtoms, dihedralAtoms, groupLabels = None):
        ''' Make structures (list) with info '''
        aseStructures = getSpeciesListCIF(structureFilename)
        
        structures, badStructures = [], []
        for indexS, s in enumerate(aseStructures):
            c = Crystal(splitASEAtomsToMols(s))
            c.calculateClusterInformation(axisAtoms, dihedralAtoms, groupLabels = groupLabels)
            #can get energypermol??
            structures.append({'crystal': c,
                               'zPrime':  len(c.asymmetricMolecules)
                               })
#            except:
#                badStructures.append({'structure': s})

        self.structures = structures
#        self.badStructures = badStructures #only include try except at very end

    def separateClusterInfo(self, _tol=1.e-4):
        outData = [[], [], [], [], [], [], [], []]
        _xVec, _yVec, _zVec = np.array([1., 0., 0.]), np.array([0., 1., 0.]), np.array([0., 0., 1.])

        for s in self.structures:
            for m in s['crystal'].asymmetricMolecules:
                for p in m.environmentPointList():
#                    print p
                    found = sum([len(x) for x in outData])
                    if abs(p['detOperation'] - 1.) < _tol:
                        if abs(p['thetaRot']) < _tol:
                            outData[0].append(p)
                    
                        elif abs(abs(np.dot(p['axisRot'], _xVec)) - 1.) < _tol:
                            outData[1].append(p)
                        elif abs(abs(np.dot(p['axisRot'], _yVec)) - 1.) < _tol:
                            outData[2].append(p)
                        elif abs(abs(np.dot(p['axisRot'], _zVec)) - 1.) < _tol:
                            outData[3].append(p)
                    else:
                        if abs(p['thetaRot']) < _tol:
                            outData[4].append(p)
                        elif abs(abs(np.dot(p['axisRot'], _xVec)) - 1.) < _tol:
                            outData[5].append(p)
                        elif abs(abs(np.dot(p['axisRot'], _yVec)) - 1.) < _tol:
                            outData[6].append(p)
                        elif abs(abs(np.dot(p['axisRot'], _zVec)) - 1.) < _tol:
                            outData[7].append(p)

        print "Total data %s and assigned %s points"%(sum([len(m.environmentPointList()) for s in self.structures for m in s['crystal'].asymmetricMolecules]), 
                                                  sum([len(x) for x in outData]))
        return outData



    
