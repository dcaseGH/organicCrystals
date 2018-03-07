import numpy as np
from ase.units import create_units
aseUnits = create_units('2014')
BOHRTOANGSTROM = aseUnits['Bohr']
RADTODEG       = 180./np.pi

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def acceptableAtomicLine(line, fileFormat='punch'):
    ''' Each atom begins with a line like this: 
        the element and a coordinate '''
    parts = line.split()
    if fileFormat == 'punch':
        if len(parts) > 3 and parts[0].isalpha() and all(map(is_number, parts[1:4])):
            return True
        else: 
            return False
    elif fileFormat == 'dma':
        # assume [number] [neighcrys string] [position] [stuff]
        if len(parts) > 4 and parts[1][0].isalpha() and all(map(is_number, parts[2:5])):
            return True
        else: 
            return False
        
class MultipoleFile():
    def __init__(self, fileName=None):
        if fileName is not None:
            self.initFromFile(fileName)
        else:
            self.aseAtoms = None
            self.atomList = None
        self.rotationInfo = None
            
#    @classmethod
    def initFromFile(self, fileName, _fileFormat = None):
        ''' Record all multipoles for these atoms 
            Assume rank 4 '''

        #formats are dma or punch
        if _fileFormat is None:
            _fileFormat = fileName.split('.')[-1].lower()

        with open(fileName, 'r') as inFile:
            lines = inFile.readlines()

        atomList = []
        for l in xrange(len(lines)):
            if _fileFormat == 'punch':
                if acceptableAtomicLine(lines[l]):
                    atomList.append({'symbol'      : lines[l].split()[0],
                                     'rank'        : 4,
                                     'position'    : np.array(map(float, lines[l].split()[1:4])),
                                     'realMoments' : [[float(lines[l+1].replace('\n', '').split()[0])],
                                                      map(float, lines[l+2].split()),
                                                      map(float, lines[l+3].replace('n', '').split()),
                                                      map(float, lines[l+4].replace('n', '').split()) +\
                                                      map(float, lines[l+5].replace('n', '').split()),
                                                      map(float, lines[l+6].replace('n', '').split()) +\
                                                      map(float, lines[l+7].replace('n', '').split())]})
            elif _fileFormat == 'dma':
                if acceptableAtomicLine(lines[l], fileFormat = _fileFormat):
                    atomList.append({'symbol'      : lines[l].split()[1][:2].split('_')[0],
                                     'rank'        : 4,
                                     'position'    :BOHRTOANGSTROM * np.array(map(float, lines[l].split()[2:5])),
                                     'realMoments' : [[float(lines[l+1].replace('\n', '').split()[0])],
                                                      map(float, lines[l+2].split()),
                                                      map(float, lines[l+3].replace('n', '').split()),
                                                      map(float, lines[l+4].replace('n', '').split()),
                                                      map(float, lines[l+5].replace('n', '').split()) +\
                                                      map(float, lines[l+6].replace('n', '').split())]})

        self.atomList = atomList

    def setAtomTypes(self, customLabels = {}):
        ''' The atom type is just the symbol (element) unless a custom label is set 
            i.e. as customLabels={22:'X'} with 22 being the index in atomList '''

        for i in xrange(len(self.atomList)):
            if i in customLabels.keys():
                self.atomList[i]['type'] = customLabels[i]
            else:
                self.atomList[i]['type'] = self.atomList[i]['symbol']
                
    def stringFormat(self, header = 'generated multipole file', rotationInfo = None, printTypes = False):
        ''' returns a string with info 
            give rotationInfo = {'axis':.., 'angle':....} if wanted - or instance can know it already
            Set the atomic Types if needed beforehand '''

        #Look to see if instance.rotationInfo has already been set
        if rotationInfo is None:
            rotationInfo = self.rotationInfo
        
        def atomString(at, rotInfo, printType):
            ''' String for punch file (N.B. change top line for .dma '''

            _outString = "%s "%at['symbol'] + " ".join(["%s"%x for x in at['position']])

            if rotInfo is not None:
                _outString += " Rotated by %s about "%(RADTODEG * rotInfo['angle']) + " ".join(["%s"%x for x in rotInfo['axis']])
            _outString += " Rank %s"%(at['rank'])

            if printType:
                _outString += " Type %s"%(at['type'])
                
            _outString += '\n'

            _outString += "%s\n"%(at['realMoments'][0][0])
            _outString += " ".join(["%s"%x for x in at['realMoments'][1]]) + '\n'
            _outString += " ".join(["%s"%x for x in at['realMoments'][2]]) + '\n'
            _outString += " ".join(["%s"%x for x in at['realMoments'][3][:5]]) + '\n'
            _outString += " ".join(["%s"%x for x in at['realMoments'][3][5:]]) + '\n'
            _outString += " ".join(["%s"%x for x in at['realMoments'][4][:5]]) + '\n'
            _outString += " ".join(["%s"%x for x in at['realMoments'][4][5:]])

            return _outString
            
        outString = '! ' + header + '\n\n'
        outString += '\n\n'.join([atomString(x, rotationInfo, printTypes) for x in self.atomList])

        return outString

    def invertRealMoments(self, mirrorPositions=True):
        ''' Uses equations from Rui's script (invertdma.f90 in orient directory on xenon) '''
        #rui inverts z, not x,y and z
        for a in self.atomList:
            if mirrorPositions:
                a['position'] *= np.array([1., 1., -1.])

            a['realMoments'][1][0] *= -1.
            a['realMoments'][2][1] *= -1.
            a['realMoments'][2][2] *= -1.
            a['realMoments'][3][0] *= -1.
            a['realMoments'][3][3] *= -1.
            a['realMoments'][3][4] *= -1.
            a['realMoments'][4][1] *= -1.
            a['realMoments'][4][2] *= -1.
            a['realMoments'][4][5] *= -1.
            a['realMoments'][4][6] *= -1.

    def copy(self, zMirror=False):
        from copy import deepcopy
        newInstance = deepcopy(self)
        if zMirror:
            newInstance.invertRealMoments()
        return newInstance

    def updatePositions(self, newPositions):
        ''' Change atomic positions to be newPositions - N.B. same order atoms '''
        for i in xrange(len(self.atomList)):
            self.atomList[i]['position'] = newPositions[i]

    def atomPositions(self):
        return np.array([x['position'] for x in self.atomList])

    def newPositions(self, newPositions, tol=1.e-3):
        ''' From base instance, give newPositions
            What is returned is a new instance with axis and angle to do rotation
            and mirrored if needed
            If not possible to overlay (with rotation) new and old positions, raise exception 
            Atomic ordering MUST BE THE SAME IN NEW POSITIONS AS OLD '''

        #Note - need safe quaternion to axis/angle
        from listMathsSlow import overlay_points_RMSD
        from quaternions   import inverse_quat_rot2Safe

        #assert that template centroid at 0. ??
        
        #get the axis and angle that relates the molecules
        originalPositions = np.array([x['position'] for x in self.atomList])
        originalPositions -= np.mean(originalPositions, axis=0)
        _newCentroid = np.mean(newPositions, axis=0)
        print 'WARNING - DOCUMENTATION IMPLIES overlay_points_RMSD OTHER WAY AROUND'
        _rmsd, _qRot = overlay_points_RMSD(originalPositions,
                                           newPositions - _newCentroid)

        #get safe inverse_quat_rot2 to handle zero angles
        if _rmsd < tol:
            angle, axis  = inverse_quat_rot2Safe(_qRot)
            newMultipole = self.copy()
        else:
            _rmsdOld = _rmsd
            # mirror not inversion
            # use a mirrored copy of old punch template
            mirrorMultipole = self.copy(zMirror=True)
            mirrorPositions = mirrorMultipole.atomPositions()
            mirrorPositions -= np.mean(mirrorPositions, axis=0)
            _rmsd, _qRot = overlay_points_RMSD(newPositions - _newCentroid,
                                               mirrorPositions)
            if _rmsd < tol:
                angle, axis = inverse_quat_rot2Safe(_qRot)
                newMultipole = self.copy(zMirror=True)
            else:
                raise Exception('Failure to overlay molecules when making new punch file - RMSD %s %s'%(_rmsdOld, _rmsd))

#        from quaternions import make_quat_rot2
#        print make_quat_rot2(angle, axis), _qRot
        newMultipole.rotationInfo = {'axis' : axis,
                                     'angle': angle}

        newMultipole.updatePositions(newPositions)
        return newMultipole

    def aseAtoms(self):
        ''' Take symbols and positions from atomList
            Function as atomList may be updated 
            Cell is just large P1 box '''
        from ase import Atoms
        from ase.spacegroup import Spacegroup
        return Atoms(symbols   = [x['symbol'] for x in self.atomList],
                     positions = [x['position'] for x in self.atomList],
                     cell      = 50. * np.eye(3),
                     info = {'spacegroup': Spacegroup(1)})
    

    def writeFile(self, filename):
        ''' Write molecule to P1 cell to check things '''
        from ase import Atoms
        from ase.spacegroup import Spacegroup
        dummyAtoms = Atoms(symbols   = [x['symbol'] for x in self.atomList],
                           positions = [x['position'] for x in self.atomList],
                           cell      = 50. * np.eye(3),
                           info = {'spacegroup': Spacegroup(1)})
        dummyAtoms.write(filename)
