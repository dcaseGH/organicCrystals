import pickle
import numpy as np

#def optimalZRotation(pts1, pts2):
#    ''' Find best overlay (which way around?) for angle around z '''
#    from scipy.optimize import fmin #clearly analytic result must be possible #powells method may be better?
#    from quaternions import quaternion_rotatn, make_quat_rot2
#    from listMathsSlow import rmsd
    

#    def functionOfOverlay(angle, pts1, pts2):
#        _quat = make_quat_rot2(angle, np.array([0., 0., 1.]))
#        rotatedPts1 = np.array([quaternion_rotatn(x, _quat) for x in pts1])
#        return rmsd(rotatedPts1, pts2)

##    optimisationOutputs = fmin(functionOfOverlay, 0., (pts1, pts2), full_output=True)
#    optimisationOutputs = fmin(functionOfOverlay, 0., (pts2, pts1), full_output=True)
##    print optimisationOutputs
#    return (make_quat_rot2(optimisationOutputs[0][0], np.array([0., 0., 1.])),
#            optimisationOutputs[1])
    
#print 'draw the molecules to check sort is picking the lowest one along c'
def calculatePerfectSurfaceAngles(surfaceCif, nanocrystalCif):
    ''' look at bottom molecule surface and return angle such that nanocrystal will match this '''

    from analyseClusters import Crystal
    from listMathsSlow   import overlay_points_RMSD
    from quaternions     import zyz_angles

    surfaceCrystal   = Crystal.fromCif(surfaceCif)
    nanoCrystal      = Crystal.fromCif(nanocrystalCif)

    #nanoCrystal is lowest along c axis
    nanoCrystal.asymmetricMolecules = [sorted(nanoCrystal.asymmetricMolecules, key = lambda x: np.dot(x.aseAtoms.get_center_of_mass(scaled=False),
                                                                                                      nanoCrystal.aseCell[2]))[0]]

    #surfaceCrystal is just lowest along c axis (dont need any more at this stage - just want to make 'perfect' cut)
    surfaceCrystal.asymmetricMolecules = [sorted(surfaceCrystal.asymmetricMolecules, key = lambda x: np.dot(x.aseAtoms.get_center_of_mass(scaled=False),
                                                                                                         surfaceCrystal.aseCell[2]))[0]]

    #should assert that both are in a unit cell 0,1
    
    #
    nanoAtoms    = nanoCrystal.asymmetricMolecules[0].aseAtoms.get_positions()[7:13,:]
    surfaceAtoms = surfaceCrystal.asymmetricMolecules[0].aseAtoms.get_positions()[7:13,:]
#    print nanoCrystal.asymmetricMolecules[0].aseAtoms.get_chemical_symbols(), surfaceCrystal.asymmetricMolecules[0].aseAtoms.get_chemical_symbols()
    print 'assuming numbering system such that ring of atoms is 2nd ring ->7-12', nanoAtoms.shape, surfaceAtoms.shape
    

    _rmsd, _qRot = overlay_points_RMSD(surfaceAtoms - np.mean(surfaceAtoms, axis=0),
                                       nanoAtoms    - np.mean(nanoAtoms, axis=0))

    assert(_rmsd < 0.1)
    return zyz_angles(_qRot)


def makeFilesSingleMolSurface(surfaceCif, nanocrystalCif, surfaceDMA, nanocrystalDMA, 
                              nanocrystalAngles = [0., 0., 0.], 
                              nanocrystalSize = [1, 1, 1],
                              freezeRotations = False, 
                              manuallyRotateNanocrystal = False,
                              rotateAllAngles = False):

    ''' Writes the files needed to run orient (.in is the input) 
        Also can write other things to visualize components or composite - may be useful as check '''

    from analyseClusters import getSpeciesListCIF, Surface, splitASEAtomsToMols, Crystal
    from orientIO import NanocrystalOnSurfaceInput
    #from ccdc.io import CrystalReader
    from ioAndInterfaces import ccdcCrystalToASE
    from multipoleFile import MultipoleFile, BOHRTOANGSTROM
    from copy import deepcopy

    #variables
#    nanocrystalAngles = [0., 0., 0.]

    #reading filenames
    moleculeFilename = 'nanoCryst.punch'

    #writing filenames
    cellFilename  = 'surface.cell'
    inputFilename = 'surfaceCalc.in' 

    #other filenams
    potentialFilename = 'fit4orient.pots'


    surfaceCrystal   = Crystal.fromCif(surfaceCif)
    nanoCrystal      = Crystal.fromCif(nanocrystalCif)

    #keep only lowest molecule along c axis
    print "%s mols in nanocrystal %s in surface"%(len(nanoCrystal.asymmetricMolecules),
                                                  len(surfaceCrystal.asymmetricMolecules))
    nanoCrystal.asymmetricMolecules = [sorted(nanoCrystal.asymmetricMolecules, key = lambda x: np.dot(x.aseAtoms.get_center_of_mass(scaled=False),
                                                                                                      nanoCrystal.aseCell[2]))[0]]
    print 'MASSIVE HACK - KEEPING BUT A SINGLE MOL'
    surfaceCrystal.asymmetricMolecules = [sorted(surfaceCrystal.asymmetricMolecules, key = lambda x: np.dot(x.aseAtoms.get_center_of_mass(scaled=False),
                                                                                                surfaceCrystal.aseCell[2]))[0]]

    # If doing this, rotate the nanocrystal so that the bottom molecule of surface is same orientation as it
    if manuallyRotateNanocrystal:

        from listMathsSlow   import overlay_points_RMSD
        from quaternions     import quaternion_rotatn

        targetAtomPositions = sorted(surfaceCrystal.asymmetricMolecules, key = lambda x: np.dot(x.aseAtoms.get_center_of_mass(scaled=False),
                                                                                                surfaceCrystal.aseCell[2]))[0].aseAtoms.get_positions()[7:13,:]

        nanoAtoms    = nanoCrystal.asymmetricMolecules[0].aseAtoms.get_positions()[7:13,:]
        print 'assuming numbering system such that ring of atoms is 2nd ring ->7-12', nanoAtoms.shape, targetAtomPositions.shape
    
        if rotateAllAngles:
            print 'dont do this - rotating all angles prevents coplanar surfaces';exit()
            _rmsd, _qRot = overlay_points_RMSD(nanoAtoms           - np.mean(nanoAtoms, axis=0),
                                               targetAtomPositions - np.mean(targetAtomPositions, axis=0))
            assert(_rmsd < 0.1)        
        else:
            #just rotate around z axis - the allowRMSDNoRestriction should not be used- just make sure atomic numbering consistent
            from listMathsSlow import optimalZRotation
            _rmsd, _qRot = optimalZRotation(nanoAtoms           - np.mean(nanoAtoms, axis=0),
                                            targetAtomPositions - np.mean(targetAtomPositions, axis=0),
                                            allowRMSDNoRestriction = False)

        def tempDistPrint(x):
            print [np.linalg.norm(x[i+1] - x[i]) for i in range(5)]

#        tempDistPrint(targetAtomPositions - np.mean(targetAtomPositions, axis=0))
#        tempDistPrint(np.array([quaternion_rotatn(x, _qRot) for x in nanoAtoms - np.mean(nanoAtoms, axis=0)]))

        #rotate nanocrystal atoms and cell
        nanoCrystal.asymmetricMolecules[0].rotateASEAtomsAroundCentroidQuat(_qRot)
        nanoCrystal.aseCell = np.array([quaternion_rotatn(x, _qRot) for x in nanoCrystal.aseCell])


    print "%s mols in nanocrystal %s in surface"%(len(nanoCrystal.asymmetricMolecules),
                                                  len(surfaceCrystal.asymmetricMolecules))

    print "Atoms in nanoCrystal moleucules " + " ".join([str(len(x.aseAtoms)) for x in nanoCrystal.asymmetricMolecules])
    print "Atoms in surface     moleucules " + " ".join([str(len(x.aseAtoms)) for x in surfaceCrystal.asymmetricMolecules])

    punchSurface     = MultipoleFile(fileName=surfaceDMA)
    punchNanocrystal = MultipoleFile(fileName=nanocrystalDMA)
    print "%s atoms in surfaceDMA"%(len(punchSurface.aseAtoms()))
    print "%s atoms in nanocrystalDMA"%(len(punchNanocrystal.aseAtoms()))
#    punchSurface.writeFile('dummyPunchSurf.cif')
#    punchNanocrystal.writeFile('dummyPunchNano.cif')

    uniqueAtomTypes = list(set(list(surfaceCrystal.uniqueElements()) +\
                               list(nanoCrystal.uniqueElements())))

    #print 'using same atoms for surface and nanocrystal'
    print 'automatically setting all H bonded to O or N to Hn type'

    inputHandler = NanocrystalOnSurfaceInput()

    #displace by c of surfaceCrystal - and add a bit
    displacementBohr = (surfaceCrystal.aseCell[2] + np.array([0., 0., 3.5])) / BOHRTOANGSTROM

    variables = [['x1',  displacementBohr[0], 'Bohr'],
                 ['y1',  displacementBohr[1], 'Bohr'],
                 ['z1',  displacementBohr[2], 'Bohr'],
                 ['alpha1', nanocrystalAngles[0], 'Degree'],
                 ['beta1',  nanocrystalAngles[1], 'Degree'],
                 ['gamma1', nanocrystalAngles[2], 'Degree']]


    atomsInCell         = [punchSurface.newPositions(m.aseAtoms.positions) for m in surfaceCrystal.asymmetricMolecules]
    atomsInNanocrystal  = [punchNanocrystal.newPositions(m.aseAtoms.positions) for m in nanoCrystal.asymmetricMolecules]

    with open(moleculeFilename, 'w') as outf: 
        #acidic hydrogens
        for im, m in enumerate(atomsInNanocrystal):
            m.setAtomTypes(dict([(x, 'Hn') for x in nanoCrystal.asymmetricMolecules[im].indicesAcidicHydrogens()]))
        outf.write('\n'.join([x.stringFormat(header = '', printTypes=True)
                              for x in atomsInNanocrystal]))
                         
    print 'passing one molecule in nanocrystal as only needs cell -- clear up later'
    aseAtomsNanocrystal = nanoCrystal.asymmetricMolecules[0].aseAtoms    

    with open(cellFilename, 'w') as outf:
        #acidic hydrogens
        for im, m in enumerate(atomsInCell):
            m.setAtomTypes(dict([(x, 'Hn') for x in surfaceCrystal.asymmetricMolecules[im].indicesAcidicHydrogens()]))
        outf.write(inputHandler.cellAndSitesString(surfaceCrystal.aseCell, atomsInCell))

    with open(inputFilename, 'w') as outf:
        outf.write(inputHandler.inputString(variables        = variables,
                                            moleculeFilename = moleculeFilename,
                                            potentialFile    = potentialFilename,
                                            cellFile         = cellFilename,
                                            nanocrystal      = aseAtomsNanocrystal,
                                            surfaceCrystal   = surfaceCrystal,
                                            nanocrystalSize  = nanocrystalSize,
                                            atomTypes        = uniqueAtomTypes,
                                            freezeRotations  = freezeRotations)
                   )

    #
    # EVERYTHING BELOW HERE IS FOR MY (DHC) DEBUGGING AND LOOKING AT FILES TO CHECK THINGS - NOT USED
    #

    #write temp.xyz with all atoms
    print 'writing temporary xyz file to temp.xyz - note that this used only by DHC, not ORIENT itself'
    try:
        from ase import Atoms
    except ImportError:
        print 'Need ASE to print out temp.xyz but not important for ORIENT - ignore this'
        return
#    superSurface = deepcopy(surfaceCrystal)

#    superSurface.filledUnitCellMolecules(superCell=np.array([3,3,2]))
#    print [x[1] for x in variables[3:6]]
    nanocrystalSize=[1,1,1];print 'hack small nano'
    surfaceBounds = inputHandler.boundsShowSurface(nanocrystalSize,
                                                   nanoCrystal.aseCell,
                                                   surfaceCrystal.aseCell,
                                                   [x[1] for x in variables[3:6]])
    print nanocrystalSize
#    print len(superSurface.asymmetricMolecules)
#    totalAtoms = Atoms(symbols = [x.symbol for m in superSurface.asymmetricMolecules for x in m.aseAtoms] +\
#                           [x.symbol for m in nanoCrystal.asymmetricMolecules for x in m.aseAtoms],
#                       positions = np.vstack([np.array([x.position for m in superSurface.asymmetricMolecules for x in m.aseAtoms]),
#                                              displacementBohr * BOHRTOANGSTROM +\
#                                                  np.array([x.position for m in nanoCrystal.asymmetricMolecules for x in m.aseAtoms])])
#                       )    

#    checkSurfPunch = MultipoleFile()
#    checkSurfPunch.initFromFile('surface.cell', _fileFormat='punch')
#    checkSurfPunch.initFromFile('dummySurface.punch', _fileFormat='punch')
#    checkNanoPunch = MultipoleFile()
#    checkNanoPunch.initFromFile('nanoCryst.punch')
#    print checkNanoPunch.atomPositions()[0]
#    print nanoCrystal.asymmetricMolecules[0].aseAtoms[0].position
#    print checkSurfPunch.atomPositions()[0]#
#    print surfaceCrystal.asymmetricMolecules[0].aseAtoms[0].position
#    tempASE = checkNanoPunch.aseAtoms()
#    tempASE.write('nanoCheck.xyz')
#    checkNanoPunch.writeFile('nanoCheck.xyz')
#    tempASE = checkNanoPunch.aseAtoms()
#    checkSurfPunch.writeFile('surfCheck.xyz')
#    exit()
#    print 'adding big gap'
    totalAtoms = Atoms(symbols = [x.symbol for m in surfaceCrystal.asymmetricMolecules for x in m.aseAtoms
                                           for a in xrange(surfaceBounds[0], surfaceBounds[1] + 1)
                                           for b in xrange(surfaceBounds[2], surfaceBounds[3] + 1)] +\
                                 [x.symbol for m in nanoCrystal.asymmetricMolecules for x in m.aseAtoms
                                           for a in xrange(0, nanocrystalSize[0])
                                           for b in xrange(0, nanocrystalSize[1])],
                       positions = np.vstack([np.array([x.position + np.dot(np.array([a, b, 1]), surfaceCrystal.aseCell)
                                                        for m in surfaceCrystal.asymmetricMolecules for x in m.aseAtoms
                                                        for a in xrange(surfaceBounds[0], surfaceBounds[1] + 1)
                                                        for b in xrange(surfaceBounds[2], surfaceBounds[3] + 1)]), 
                                              displacementBohr * BOHRTOANGSTROM +\
#                                              np.array([0., 0., 50.]) + displacementBohr * BOHRTOANGSTROM +\
                                              np.array([x.position + np.dot(np.array([a, b, 1]), nanoCrystal.aseCell)
                                                        for m in nanoCrystal.asymmetricMolecules for x in m.aseAtoms
                                                        for a in xrange(0, nanocrystalSize[0])
                                                        for b in xrange(0, nanocrystalSize[1])])])
                                                         
                       )
#    print np.min(np.array([x.position + np.dot(np.array([a, b, 1]), surfaceCrystal.aseCell)
#                                                        for m in surfaceCrystal.asymmetricMolecules for x in m.aseAtoms
#                                                        for a in xrange(surfaceBounds[0], surfaceBounds[1] + 1)
#                                                        for b in xrange(surfaceBounds[2], surfaceBounds[3] + 1)]), axis=0)
#    print np.max(np.array([x.position + np.dot(np.array([a, b, 1]), surfaceCrystal.aseCell)
#                                                        for m in surfaceCrystal.asymmetricMolecules for x in m.aseAtoms
#                                                        for a in xrange(surfaceBounds[0], surfaceBounds[1] + 1)
#                                                        for b in xrange(surfaceBounds[2], surfaceBounds[3] + 1)]), axis=0)
#    print np.min(                                              displacementBohr * BOHRTOANGSTROM +\
#                                              np.array([x.position + np.dot(np.array([a, b, 1]), nanoCrystal.aseCell)
#                                                        for m in nanoCrystal.asymmetricMolecules for x in m.aseAtoms
#                                                        for a in xrange(0, nanocrystalSize[0])
#                                                        for b in xrange(0, nanocrystalSize[1])]), axis=0)
#    print np.max(                                              displacementBohr * BOHRTOANGSTROM +\
#                                              np.array([x.position + np.dot(np.array([a, b, 1]), nanoCrystal.aseCell)
#                                                        for m in nanoCrystal.asymmetricMolecules for x in m.aseAtoms
#                                                        for a in xrange(0, nanocrystalSize[0])
#                                                        for b in xrange(0, nanocrystalSize[1])]), axis=0)
#    print displacementBohr * BOHRTOANGSTROM
#    print surfaceCrystal.aseCell[2]
    totalAtoms.write('temp.xyz')
