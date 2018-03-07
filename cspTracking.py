''' Objects to help track CSP structures through CrystalPredictor and CrystalOptimiser '''
from ase.units import create_units
aseUnits = create_units('2014')
HARTREE2KJPM = (aseUnits['Hartree']*aseUnits['mol'])/aseUnits['kJ']

def generalCallErrorCatch(instance, function, argument):
    ''' Allows instance.function(argument) with generic errors '''
    try:
        return getattr(instance, function)(argument)
    except:
        return None

def processFile(summaryFilePointer, inDict):
    ''' read file, inDict has the info to look for particular lines and
        processes info (see below) '''

    import os
    
    outDict = {}

    if not os.path.isfile(summaryFilePointer):
        return dict([[k, None] for k in inDict.keys()])

    #for some reason the order of loops was wrong originally
    #also can only run readlines once, so save it
    with open(summaryFilePointer, 'r') as inFile:
        savedList = inFile.readlines()
        for k in inDict.keys():
            try:
                for l in savedList:
                    if inDict[k][0] in l:
                        outDict[k] = inDict[k][2]((l.split()[inDict[k][1]]))
            except:
                outDict[k] = None

    # add entry = None for things that weren't picked up
    for k in inDict.keys():
        if k not in outDict.keys():
            outDict[k] = None
            
    return outDict

#possibly too general, move if need
def processSummary(summaryFilePointer, inDict, calculateZp=True):
    ''' This reads a file and returns a dictionary with data as per instructions in inDict
        inDict is {'common key for in and out dict': [string to look for, index of number after split, data type]} 
        returns outDict{common key for in and out dict: piece of data} '''

    import os
    
    outDict = {}

    if not os.path.isfile(summaryFilePointer):
        return dict([[k, None] for k in inDict.keys()])

    with open(summaryFilePointer, 'r') as inFile:
        for l in inFile.readlines():
            for k in inDict.keys():
                try:
                    if inDict[k][0] in l:
                        outDict[k] = inDict[k][2]((l.split()[inDict[k][1]]))
                except:
                    outDict[k] = None

    #calculate z' (all whole molecules) and one day intra energy if needed
    if calculateZp:
        with open(summaryFilePointer, 'r') as inFile:
            summaryString = inFile.read().split('Molecule -- Intramolecular Energy contribution (kJ/mol)')[1]\
                                         .split('Total lattice energy (kJ/mol)')[0]
            outDict['zPrime'] = len([x for x in summaryString if x=='\n']) - 1
    return outDict


class CSPStructure():
    ''' Object to make it easier to manage structures '''

    def __init__(self,
                 idNumber,
                 initCrystalOptimiserData = None,
                 initCrystalPredictorData = None,
                 initPCMData              = None):

        ''' Initialize with an idNumber - possibly also pass a pointer to files from which
            CrystalOptimiser or CrystalPredictor data can be taken '''
        
        self.idNumber = idNumber

        if initCrystalPredictorData:
            if type(initCrystalPredictorData) == dict:
                self.crystalPredictorData = initCrystalPredictorData
            elif type(initCrystalPredictorData) == str:
                self.createCrystalPredictorData(initCrystalOptimiserData)
            else:
                raise TypeError("initCrystalPredictorData must be a dict of data or str pointer")
        else:
            self.crystalPredictorData = {} 

        if initCrystalOptimiserData:
            self.createCrystalOptimiserData(initCrystalOptimiserData)
        else:
            self.crystalOptimiserComplete = None #for now just mark this as None
            self.crystalOptimiserData   = {} 

        if initPCMData:
            if 'gaussianLog' in initPCMData.keys():
                self.pcmData = self.createDMACRYSData(initPCMData['runDirectory'],
                                                      gaussianLog = initPCMData['gaussianLog'])
            else:
                self.pcmData = self.createDMACRYSData(initPCMData['runDirectory'])
        else:
            self.pcmData = {}

    def createDMACRYSData(self,
                          runDirectory,
                          finalStructureFile = 'fort.16',
                          summaryFile        = 'fort.12',
                          calculation        = 'pcm',
                          gaussianLog        = None,
                          templateDict       = {'eInterFinal': ['Final Lattice Energy:', 3, float]}):

        ''' Assumes DMACRYS run in a directory - and if gaussianLog not None, tries to get
            energy from here '''

        import os

        currentData = processSummary(runDirectory + '/' + summaryFile,
                                     templateDict,
                                     calculateZp = False)

        currentData['finalStructurePointer'] = os.path.abspath(runDirectory + '/' + finalStructureFile)
        currentData['finalStructureFileExists'] =\
                    os.path.isfile(currentData['finalStructurePointer'])

        #eTotalFinal is per mol, but components (intra and inter) are not (for now dhc 251017)
        if gaussianLog is not None:
            if '*' in gaussianLog: #could always just glob a single file?? save code?
                import glob
                for imol, gl in enumerate(sorted(glob.glob(gaussianLog))):
                    currentData.update(processSummary(gl,
                                                      {'eIntraHartrees_mol%s'%(imol+1): ['SCF Done:', 4, float]},
                                                      calculateZp = False))

                currentData.update({'eIntraHartrees':
                                    sum([currentData[x] for x in currentData.keys() if 'eIntraHartrees_mol' in x])})
                _tempZPrime = len(glob.glob(gaussianLog))
            else:
                currentData.update(processSummary(gaussianLog,
                                                  {'eIntraHartrees': ['SCF Done:', 4, float]},
                                                  calculateZp = False))
                _tempZPrime = 1
                
            try:
                currentData['eIntraFinal'] = currentData['eIntraHartrees'] * HARTREE2KJPM
                currentData['eTotalFinal'] = (currentData['eInterFinal'] + currentData['eIntraFinal']) / float(_tempZPrime)
#                print currentData, _tempZPrime, 
            except:
                currentData['eIntraFinal'], currentData['eTotalFinal'] = None, None

                
        if calculation == 'pcm':
            try: 
                # pcmData = currentData is needed for packingCoefficient to point to structure
                self.pcmData = currentData
                currentData['packingCoefficient'] = self.packingCoefficient(stage='pcmData')
            except (RuntimeError, IOError) as error:
                print 'ccdc failed for ', self.idNumber, error

        return currentData

    def createCrystalOptimiserData(self,
                                   crystalOptimiserSummaryFile,
                                   templateDict = {'eInterFinal': ['Intermolecular energy contribution (kJ/mol)', 4, float],
                                                   'eTotalFinal': ['Total lattice energy (kJ/mol)',               4, float],
                                                   'rhoFinal':    ['Crystal density (g/cm3)',                     3, float]},
                                   lowestStructurePointer = 'lowest.res'
                                   ):

        ''' At the moment, this just gets standard data from summary.out 
            change templateDict if you want different information
            general format is key: [string to look for, column data is in, type]
            If calculation has failed/not been run - can still point to where summary.out would be '''

        import os

        #at moment pass a dummy summary.out in if doesnt exist,, in which case dont save

        if os.path.isfile(crystalOptimiserSummaryFile):
            self.crystalOptimiserSummaryFile = crystalOptimiserSummaryFile
        else:
            self.crystalOptimiserSummaryFile = None

        #this is not the best - lowest structure pointer should be crystalOptimiserData['lowestStructurePointer']
        self.lowestStructurePointer = lowestStructurePointer
        self.crystalOptimiserDir  = '/'.join(os.path.realpath(crystalOptimiserSummaryFile).split('/')[:-1])
        self.crystalOptimiserData = processSummary(crystalOptimiserSummaryFile,
                                                   templateDict)
        self.crystalOptimiserData['runDirectory'] = '/'.join(os.path.realpath(crystalOptimiserSummaryFile).split('/')[:-1])
        
        self.crystalOptimiserComplete, self.crystalOptimiserError = self.crystalOptimiserLogFinished()
#        self.crystalOptimiserData['lowestStructurePointer'] = lowestStructurePointer
#        self.crystalOptimiserData['lowestStructureFileExists'] =\
#             os.path.isfile(self.crystalOptimiserDir + '/' + lowestStructurePointer)
        self.crystalOptimiserData['finalStructurePointer'] = self.crystalOptimiserData['runDirectory'] + '/' + lowestStructurePointer
        self.crystalOptimiserData['finalStructureFileExists'] =\
             os.path.isfile(self.crystalOptimiserData['finalStructurePointer'])
        if self.crystalOptimiserData['finalStructureFileExists']:
            self.crystalOptimiserData['packingCoefficient'] = self.packingCoefficient(stage='crystalOptimiserData')
        
        return #error check??

    def crystalOptimiserLogFinished(self, logName = '/CrystOpt.log'):
        ''' Reads a log and checks last lines to see if it worked '''

        import os
        if not os.path.isfile(self.crystalOptimiserDir + logName):
            return (False, 'Cant find file')
        
        with open(self.crystalOptimiserDir + logName, 'r') as sumFile:
            lastLines = "".join(sumFile.readlines()[-2:])
            if 'CrystalOptimizer terminated successfully' in lastLines:
                return (True, None)
            elif 'Check the compatibility of the electrostatic LAM database with internal record. Exiting.' in lastLines:
                return (False, "electroStaticDBError_%s"%lastLines.split('in the electrostatic LAM database')[1]\
                                                                  .split('Check the compatibility')[0].strip().replace('\n',''))
            elif 'Check the compatibility of the intramolecular database with internal record. Exiting.' in lastLines:
                return(False, "intramolecularDBError_%s"%lastLines.split('in the intramolecular database entry')[1]\
                                                                 .split('Check the compatibility')[0].strip().replace('\n', ''))
        return (False, 'unknown error')

    def readFreeEnergyData(self,
                           freeEnergyLog,
                           infoDict = {'ZPE' : ['Neat zero point energy',  5, float],
                                       'Evib': ['Neat vibrational energy', 4, float]}):
        ''' freeEnergyLog is a string pointer to the log file
            extracts information that may be useful 
            Evib is the Helmholtz energy per molecule (inc. ZPE) from harmonic approx,
            so without debye, and without smearing type factors (KDE) '''

        import os
        
        freeEnergyData                  = processFile(freeEnergyLog, infoDict)
        freeEnergyData['logPointer']    = os.path.abspath(freeEnergyLog)
        freeEnergyData['logFileExists'] = os.path.isfile(freeEnergyLog)

        self.freeEnergyData = freeEnergyData
        
    def setAttribute(self, attribute):
        ''' This is like getAttribute except pass in Attribute instance with pointer and value 
            Have pointer in same format as getAttribute '''

        at, newValue = attribute.pointer, attribute.value
        if type(at) == str:
            # property attributes dont appear in __dict__ so use getattr (could always use this??)
            if at not in self.__dict__.keys():
                try:
                    setattr(self, at, newValue)
                    return getattr(self, at)
                except ValueError:
                    return None
            else:
                self.__dict__[at] = newValue
                return self.__dict__[at]

        if type(at) == dict:
            if type(at[at.keys()[0]]) == str:
                self.__dict__[at.keys()[0]][at[at.keys()[0]]] = newValue
                return self.__dict__[at.keys()[0]][at[at.keys()[0]]]
        return None
             
    def getAttribute(self, at):
        ''' If at is a string, return self[at], if a {key:value}, self[key[value]]
            nesting up to two dictionaries is possible 
            This works for self.attributes which are dict ie crystalOptimiserData etc '''

        if type(at) == str:
            # property attributes dont appear in __dict__ so use getattr (could always use this??)
            if at not in self.__dict__.keys():
                try:
                    return getattr(self, at)
                except ValueError:
                    return None
                except AttributeError:
                    return None
            else:
                return self.__dict__[at]

        if type(at) == dict:
            if type(at[at.keys()[0]]) == str:
                #this kind of line would be if  self.__dict__[at.keys()[0]] was not a dict but a class
#                return self.__dict__[at.keys()[0]].__dict__[at[at.keys()[0]]]
                return self.__dict__[at.keys()[0]][at[at.keys()[0]]]
#                try:
#                    return self.__dict__[at.keys()[0]][at[at.keys()[0]]]
#                except KeyError:
#                    return None
#       do this another time if needed     elif type(at[at.keys()[0]]) == dict:
#                if type(at[at.keys()[0]]) == str:
#                    return at[at.keys()[0]]
        return None
   #certain properties may need @property - 

    def meetsConditions(self, requiredAttributes = []):
        ''' First check self has the attributes demanded, then possibly further checks on limits '''

        for a in requiredAttributes:

            try:
                tempValue = self.getAttribute(a.pointer)
                if a.value is not None:
                    assert(tempValue == a.value)
                if a.notValue is not None:
                    assert(tempValue != a.notValue)
                if a.minValue is not None:
                    assert(tempValue > a.minValue)
                if a.maxValue is not None:
                    #in Python None is always less than a value (even a negative one)
                    assert(tempValue is not None)
                    assert(tempValue < a.maxValue)
            except:
                return False
        return True

#    def toCCDCStructure(self, pointer = None, stage = 'crystalOptimiserData', identifier = None):
    def toCCDCStructure(self, pointer = None, stage = None, identifier = None):
        ''' Returns CCDC structure instance from appropriate pointer
            default atm is crystal optimiser 
            stage is SetCSPStructures.calculation = 'Data' (as general rule) 
            give identifier to change the name in the cif '''

        from ccdc.io import CrystalReader

        # temporary fail as code work in progress
        if stage not in ['crystalOptimiserData', 'pcmData']:
            raise Exception('Implement other options to get structures if needed')

        #generalize this to be self.getAttribute(pointer = {stage: lowestStructurePointer}) if pointer is None
        if pointer is None:
            pointer = self.getAttribute({stage: 'finalStructurePointer'})

        # CrystalReader can be told explicitly which file format to look for
        if pointer[-3:].lower() in ['res', '.16']: #.16 is fort.16 from DMACRYS
            fileFormat = 'res'
        elif pointer[-3:].lower() == 'cif':
            fileFormat = 'cif'

        ccdcCrystal = CrystalReader(pointer, format = fileFormat)[0]

        if identifier:
            ccdcCrystal.identifier = identifier

        return ccdcCrystal

    @property
    def cryOptEnergyPerMol(self):
        try:
            return self.getAttribute({'crystalOptimiserData': 'eTotalFinal'}) /\
                   float(self.getAttribute({'crystalOptimiserData': 'zPrime'}))
        except:
            raise ValueError('Need to set crystalOptimiserData')

    @property
    def pcmEnergyPerMol(self):
        ''' Einter + Eintra / z' (assuming have set all these things) '''
        try:
            return (self.getAttribute({'pcmData': 'eInterFinal'}) +\
                    self.getAttribute({'pcmData': 'eIntraFinal'})) /\
                    float(self.getAttribute({'pcmData': 'zPrime'}))
        except:
            raise ValueError('Need to set pcmData')
       
#    @property #technically this is a property of a property (as CCDC.crystal has property of this too??)
#    def packingCoefficient(self, stage = 'crystalOptimiserData'): #DO NOT PUT MUTABLE OBJECTS AS DEFAULT VALUE!!!!!
    def packingCoefficient(self, stage = None):
        ''' Generate CCDC structure with (crystal optimiser) data and use their algorithm
            How does it define volume of a molecule??? '''

        ccdcStructure = self.toCCDCStructure(stage=stage)
        return ccdcStructure.packing_coefficient

    def sameStructureAs(self, anotherStructure, stage=None, nRMSDMatches = 15, rmsdTol = 0.3, allowMolecularDifferences=False):
        ''' Takes final structures from a particular stage (e.g. pcmData) and 
            uses CCDC overlay to match structures within tolerances - returns True/False '''

        from ccdcPackingSimilarityTools import matchStructureLists

        matchData = matchStructureLists([self.toCCDCStructure(stage=stage)],
                                        [anotherStructure.toCCDCStructure(stage=stage)],
                                        allowMolecularDifferences=allowMolecularDifferences,
                                        returnRMSD   = True,
                                        nRMSDMatches = nRMSDMatches)

        return matchData[0][0,0] == nRMSDMatches and matchData[1][0,0] <= rmsdTol
        
class SetCSPStructures():
    ''' List of CSPStructures and subroutines to manipulate them 
        calculation will mean that structures can be from experiment/crystalPredictor/crystalOptimiser
        etc (use these exact names) and will have a default to point to for eg res file or properties '''

    def __init__(self, structures = [], calculation = None):
        self.structures  = structures

        if calculation not in [None, 'experiment', 'crystalPredictor', 'crystalOptimiser', 'pcm']:
            raise ValueError('Incorrect type or spelling of calculation in SetCSPStructures')
        self.calculation = calculation
        
    @classmethod
    def initFromCSV(cls,
                    csvFile,
                    calculations          = 'crystalOptimiser',
                    selectIDNumberLambda  = None,
                    idNumberApplyLambda   = None,
                    idNumberIndex         = None,
                    finalEnergyIndex      = None,
                    finalInterEnergyIndex = None,
                    finalDensityIndex     = None,
                    finalStructureDir     = None,
                    finalStructurePointer = None,
                    assertZPrime          = None,
                    useHeader             = False,
                    delim_whitespace      = True,
                    specifiedStudyID      = None):

        ''' Read comma separated value file (whitespace as delim, comma if delim_whitespace = False) 
            assume that all finalStructures are named [idNumber].res
            Change things around if calculations!=crys opt '''

        import numpy as np
        import pandas as pd
        import os
        
        #if the table has a header, use it unless other instructions given (to assign columns)
        if useHeader:
            dataTable = pd.read_csv(csvFile,
                                    delim_whitespace = delim_whitespace)
            columns   = [c for c in dataTable.columns]
            # These are standard labels which are used in the online data
            # note - cannot be too clever here as some labels are substrings of others
            idNumberIndex    = 'Label'

            if 'Intermolecular Lattice Energy /kJmol-1' in columns:
                finalInterEnergyIndex = 'Intermolecular Lattice Energy /kJmol-1'
            finalEnergyIndex = 'Lattice Energy /kJmol-1'
            studyIDIndex     = 'Study_ID'

            #if you are doing things from a database table ensure that you arent mixing studies
            if len(set(dataTable['Study_ID'])) > 1:
                if specifiedStudyID:
                    dataTable = dataTable.loc[dataTable['Study_ID'] == specifiedStudyID]
                else:
                    raise NameError('Specify Study_ID in initFromCSV')
        else:
            dataTable = pd.read_csv(csvFile,
                                    header           = None,
                                    delim_whitespace = delim_whitespace)

        if calculations not in ['crystalOptimiser', 'pcm']:
            raise Exception('implement other methods later')

        dataPointer = calculations + 'Data'

        # can provide a function to select data by idNumber (e.g. dont take ones with E in)
        if callable(selectIDNumberLambda):
            dataTable = dataTable.loc[dataTable['Label'].apply(selectIDNumberLambda)]
            
        #note, must use 'is not None' to catch cases where = 0
        if idNumberIndex is not None:
            idNumbers = dataTable[idNumberIndex]
        else:
            idNumbers = np.array(range(1, dataTable.shape[0]+1))

        #make a set of structures
        newInstance = cls([CSPStructure(i) for i in idNumbers])

        # change the idNumbers to be ints (so 'dfA123' -> 123 for e.g.)
        if idNumberApplyLambda:
            newInstance.applyLambdaFunctionIDNumbers(idNumberApplyLambda)
        
        #get various properties if available
        if finalEnergyIndex is not None:
            for ie, e in enumerate(dataTable[finalEnergyIndex]):
                newInstance.structures[ie].setAttribute(Attribute(pointer = {dataPointer: 'eTotalFinal'},
                                                                  value   = e))

        if finalInterEnergyIndex is not None:
            for ie, e in enumerate(dataTable[finalInterEnergyIndex]):
                newInstance.structures[ie].setAttribute(Attribute(pointer = {dataPointer: 'eInterFinal'},
                                                                  value   = e))

        if finalDensityIndex is not None:
            for ie, e in enumerate(dataTable[finalDensityIndex]):
                newInstance.structures[ie].setAttribute(Attribute(pointer = {dataPointer: 'rhoFinal'},
                                                                  value   = e))
                
        #use this if all the structures are in [directory]/[idNumber].res
        if finalStructureDir is str:
            for s in newInstance.structures:
                s.crystalOptimiserDir    = finalStructureDir
                s.lowestStructurePointer = str(s.idNumber) + '.res'

        elif callable(finalStructureDir): #i.e. a lambda function (of the idNumber)
            for s in newInstance.structures:
                s.setAttribute(Attribute(pointer = {dataPointer: 'runDirectory'},
                                         value   = finalStructureDir(s.idNumber)))
        # also set structure pointer
        if callable(finalStructurePointer):
            for s in newInstance.structures:
                s.setAttribute(Attribute(pointer = {dataPointer: 'finalStructurePointer'},
                                         value   = finalStructurePointer(s.idNumber)))
                s.setAttribute(Attribute(pointer = {dataPointer: 'finalStructureFileExists'},
                                         value   = os.path.isfile(s.getAttribute({dataPointer: 'finalStructurePointer'}))))
                try:
                    s.setAttribute(Attribute(pointer = {dataPointer: 'packingCoefficient'},
                                             value   = s.packingCoefficient(stage=dataPointer)))
                except (RuntimeError, IOError) as error:
                    print 'ccdc failed for ', s.idNumber, error
        if assertZPrime is not None:
            newInstance.setGlobalAttributes([Attribute(pointer = {dataPointer: 'zPrime'},
                                                       value   = assertZPrime)])

        newInstance.calculation = calculations
            
        return newInstance
    
    @classmethod
    def initFromCrystalPredictorLog(cls,
                                    crystalPredictorLog,
                                    crystalOptimiserSummaryPointer = None,
                                    pcmSummaryPointer              = None):

        ''' crystalOptimiserSummarypointer is a glob friendly thing with which to 
            add CO information -- assume form [base dirs]/x-1/[name of summaryfile]
            to point to thing which will have idNumber x '''
        
        from crystalPredictorSubroutines import analyseEditedLog

        newInstance = cls([CSPStructure(ix + 1, initCrystalPredictorData = x)
                           for ix, x in enumerate(analyseEditedLog(crystalPredictorLog))])

        #also add all the crystalOptimiserData/PCM if requested/possible
        if crystalOptimiserSummaryPointer:
            newInstance.addInfoFromGlobString(crystalOptimiserSummaryPointer,
                                              calculation = 'crystalOptimiser')
        if pcmSummaryPointer:
            newInstance.addInfoFromGlobString(pcmSummaryPointer,
                                              calculation = 'pcm')

        return newInstance
    
#    def addCrystalOptimiserInfoFromGlobString(self, globPointer, dummySummaryFile = 'summary.out'): #old name
    def addInfoFromGlobString(self, globPointer, dummySummaryFile = None, calculation = 'crystalOptimiser'):
        ''' Assumes that either glob pointer is of form */x-1/summaryName 
            or dirs in which case it will pretend that there is a summary in each dir 
            the position of idNumber x is important (x-1 is idNumber x s dir) '''

        import glob
        import os

        if dummySummaryFile is None:
            dummySummaryFile = {'pcm': 'fort.12',
                                'crystalOptimiser': 'summary.out'}[calculation]
        
        listOfPointers = glob.glob(globPointer)

        if any([os.path.isdir(x) for x in listOfPointers]) and not\
           all([os.path.isdir(x) for x in listOfPointers]):
            raise Exception('point purely at directories to add crystal optimiser info')
           
        # if pointing at dirs, add a dummy summary name
        if all([os.path.isdir(x) for x in listOfPointers]):
            listOfPointers = map(lambda x: x + '/' + dummySummaryFile, listOfPointers)

        for s in listOfPointers:
            currentStructure = [x for x in self.structures\
                                if x.idNumber == int(s.split('/')[-2].split('-')[0])][0]
            if calculation == 'crystalOptimiser':
                currentStructure.createCrystalOptimiserData(s)
            elif calculation == 'pcm':
                # assuming z'=1 and mol1/act_mult1.log is the file for the intra energy
                currentStructure.createDMACRYSData(s.replace('/'+dummySummaryFile, ''),
                                                   gaussianLog = s.replace(dummySummaryFile,
                                                                           'mol1/act_mult1.log'))
            else:
                raise ValueError('addInfoFromGlobString only takes pcm or crystalOptimiser as calculation arguments')

    def addFreeEnergyData(self, lambdaFunctionToLogFiles):
        ''' Give a lambda function (of idNumber) to the log files (possibly better than glob??) 
            e.g. lambda x: [myDirectory]/str(x)-1/free_energy.log '''

        # have a subset, presumably (e.g. idNumber < 2000)
        for s in self.structures:
            s.readFreeEnergyData(lambdaFunctionToLogFiles(s.idNumber))

    def toPandasDataFrame(self, dataRequirements = [], allowMissingData = False):
        ''' For now just has standard column headers (its a dict with key being name of header 
            and value being how to get it '''

        import pandas as pd

        #can have a standard list of requirements for laziness
        if len(dataRequirements) == 0:
            
            dataRequirements = [Attribute(name    = 'ID',
                                          pointer = 'idNumber'),
                                Attribute(name    = 'eTotalFinal',
                                          pointer = {'crystalOptimiserData': 'eTotalFinal'}),
                                Attribute(name    = 'rhoFinal',
                                          pointer = {'crystalOptimiserData': 'rhoFinal'}),
                                Attribute(name    = 'SG',
                                          pointer = {'crystalPredictorData': 'spaceGroup'})]


        #presumably have a max no of points
        if allowMissingData:
            return pd.DataFrame(dict([[a.name, [generalCallErrorCatch(x, 'getAttribute', a.pointer) for x in self.structures]]
                                      for a in dataRequirements]))
        else:
            #can subset to only have energies etc
            subsetStructures = self.subsetStructures(dataRequirements)

            return pd.DataFrame(dict([[a.name, [x.getAttribute(a.pointer) for x in subsetStructures]] for a in dataRequirements]))

    def applyLambdaFunctionIDNumbers(self, lambdaFunction):
        ''' Apply a lambda function to the idNumber (could extend to all properties quite easily) '''

        for i in xrange(len(self.structures)):
            self.structures[i].idNumber = lambdaFunction(self.structures[i].idNumber)
        
    # should be called subset list or something
    def subsetStructures(self, requiredAttributes = [], sortAttribute = None):
        ''' Returns a list with a subset of the structures
            Consider memory issues (pointers vs new objects) - test these things thoroughly '''

        subset = [x for x in self.structures if
                  x.meetsConditions(requiredAttributes=requiredAttributes)]

        if sortAttribute is not None:
            return sorted(subset, key=lambda x: x.getAttribute(sortAttribute.pointer))
        else:
            return subset

    def structureIDNumber(self, idNumber):
        ''' easy return of a particular idNumber 
            could always have an assert to check unique ID '''
        return self.subsetStructures(requiredAttributes=([Attribute(pointer = 'idNumber',
                                                                   value   = idNumber)]))[0] 
        
    def matchToStructureList(self,
                             listCCDCStructures,
                             stage = None,
                             requiredAttributes = [], #USE IMMUTABLE INPUTS!!!
                             returnRMSD=False,
                             allowMolecularDifferences=False):

        ''' Match current structures to this list
            returns a dict which should be saved and manipulated
            CAVEAT - potential problems with Python pointers; deepcopy if needed to use in complicated manner '''

        from ccdcPackingSimilarityTools import matchStructureLists

        #can restrict search, e.g. low energy structures
        subsetStructures = self.subsetStructures(requiredAttributes)
        data = matchStructureLists(listCCDCStructures,
                                   [x.toCCDCStructure(stage=stage) for x in subsetStructures],
                                   returnRMSD=returnRMSD,
                                   allowMolecularDifferences=allowMolecularDifferences)
        if returnRMSD:
            matchesMatrix, rmsdMatrix = data[0], data[1]
        else:
            matchesMatrix, rmsdMatrix = data, None
        
        return {'targetList'   : listCCDCStructures,
                'candidateList': subsetStructures,
                'matchesMatrix': matchesMatrix,
                'rmsdMatrix'   : rmsdMatrix}

    def calculatePackingCoefficients(self, stage=None):
        ''' Calculate pack coeff for a particular stage '''

#        if stage is None and self.calculations is not None:
#            stage = self.calculations

        #perhaps consider having a subset, as some will not have valid res files to point to?
        for s in self.structures:
            try:
                s.setAttribute(Attribute(pointer = {stage: 'packingCoefficient'},
                                         value   =  s.packingCoefficient(stage=stage)))
            except:
                print s.idNumber, 'failed'
        
    def minValueInSet(self, attributeOfInterest, safeMode=False):
        ''' Pass in an Attribute with a pointer - perhaps have maxValue or conditions if needed later 
            beware None < any number in python 
            safeMode allows things like dict.keys that haven't been defined '''

        if safeMode:
            values = [generalCallErrorCatch(x, 'getAttribute', attributeOfInterest.pointer) for x in self.structures]
        else:
            values = [x.getAttribute(attributeOfInterest.pointer) for x in self.structures]

        return min([x for x in values if x is not None])

    def setGlobalAttributes(self, listAttributesOfInterest):
        ''' Set all structures (e.g. Attribute(pointer = {'crystalOptimiserData': 'zPrime'},
                                               value   = 1) '''

        #must give a list --- why does this not work (setting non list to list I mean)???
#        print listAttributesOfInterest
#        if type(listAttributesOfInterest) is not list:
#            listAttributesOfInterest = list(listAttributesOfInterest)

        for a in listAttributesOfInterest:
            for s in self.structures:
                s.setAttribute(a)

    def statusString(self, calculation = None):
        ''' Just return a string containing numbers of successes and errors '''

        # Only look at data from a particular calculation
        if calculation is None:
            calculation = self.calculation

        # build up a list of things to say
        dataList = []
        dataList.append("Total number of all structures       %s"%(len(self.structures)))
        if calculation == 'crystalOptimiser':
            dataList.append('Crystal Optimiser Data')
            #currently do not give notValue=None
#            dataList.append("Total number of directories assigned %s"%(len(self.subsetStructures(
#                                                                       requiredAttributes=[Attribute(pointer='crystalOptimiserDir',
#                                                                                                     notValue=None)]))))
            dataList.append("Total number of lowest.res           %s"%(len(self.subsetStructures(
                            requiredAttributes = [Attribute(pointer={'crystalOptimiserData': 'finalStructureFileExists'},
                                                            value=True)]))))
            dataList.append("Total number of jobs successful      %s"%(len([x for x in self.structures if x.crystalOptimiserComplete])))
            dataList.append("Error messages " + " ".join(map(str, set([x.getAttribute('crystalOptimiserError')
                                                                       for x in self.structures]))))            
            
        return '\n'.join(dataList)

    def clusterSet(self, stage = None, nRMSDMatches = 15, rmsdTol = 0.3, clusterEnergyLimit = None, clusterEnergyTolerance = 5.):
        ''' default tolerance 15 mol overlay with 0.3 angstrom rmsd for CCDC structure overlay
            stage = 'pcmData' for e.g.
            attaches to each structure a pointer to a lower energy one that matches - these can be filtered later '''
        import numpy as np
        
        assert(stage in ['pcmData', 'crystalOptimiserData'])
        #include cutoff in subset definition
        orderedStructures = self.subsetStructures(requiredAttributes = [Attribute(pointer={stage: 'finalStructureFileExists'},
                                                                                value=True)],
                                                  sortAttribute      = Attribute(pointer={stage: 'eTotalFinal'}))

        #define a cutoff beyond which we don't cluster
        if clusterEnergyLimit is not None:
            _clusterCutoff = orderedStructures[0].getAttribute({stage: 'eTotalFinal'}) + clusterEnergyLimit
        else:
            _clusterCutoff = np.inf

        #loop through template (lower in energy) and potential matches
        for iTemIndx, templateStructure in enumerate(orderedStructures):

            #only bother with low energy structures
            if templateStructure.getAttribute({stage: 'eTotalFinal'}) > _clusterCutoff:
                break

            #has this structure already been matched to something? if so skip it
            if generalCallErrorCatch(templateStructure, 'getAttribute', {stage: 'clusterMatch'}) is not None:
                continue

            for targetStructure in orderedStructures[(iTemIndx+1):]:

                if targetStructure.getAttribute({stage: 'eTotalFinal'}) > _clusterCutoff:
                    break

                #use generalCallErrorCatch as clusterMatch may not exist
                if generalCallErrorCatch(targetStructure, 'getAttribute', {stage: 'clusterMatch'}) is not None:
                    continue

                if targetStructure.getAttribute({stage: 'eTotalFinal'}) >\
                templateStructure.getAttribute({stage: 'eTotalFinal'}) + clusterEnergyTolerance:
                    break

                # for now limit to one match (ignore possibility of matching multiple structures)
                if templateStructure.sameStructureAs(targetStructure,
                                                     stage=stage,
                                                     nRMSDMatches = nRMSDMatches,
                                                     rmsdTol = rmsdTol,
                                                     allowMolecularDifferences=False):

                    targetStructure.setAttribute(Attribute(pointer = {stage: 'clusterMatch'},
                                                           value   = templateStructure))

            templateStructure.setAttribute(Attribute(pointer = {stage: 'clusterMatch'},
                                                     value   = None))

        return # something??

    def makeCIF(self, fileName, requiredAttributes = [], sortAttribute = None):
        ''' writes a file, called fileName, with the subset of structures 
            defined by requiredAttributes 
            Make sure you set calculation so that toCCDCStructure knows which structure to take '''

        from ioAndInterfaces import ccdcCrystalsToCIF

        #change the identifier tag for things which dont have a final energy
        assert(self.calculation in ['pcm', 'crystalOptimiser'])

        #sort structures too ??? use energy in identifier if possible (round to 5 figures)
        if all([x.meetsConditions(requiredAttributes = [Attribute(pointer = {self.calculation + 'Data': 'eTotalFinal'})])\
                for x in self.subsetStructures(requiredAttributes)]):
            ccdcCrystalsToCIF([x.toCCDCStructure(stage=self.calculation + 'Data',
                                                 identifier = "%s_%s"%(x.idNumber,
                                                                       float("%.5g" % x.getAttribute(\
                                                                            {self.calculation + 'Data': 'eTotalFinal'}))))
                               for x in self.subsetStructures(requiredAttributes, sortAttribute=sortAttribute)],
                              fileName)
        else:
            ccdcCrystalsToCIF([x.toCCDCStructure(stage=self.calculation + 'Data',
                                                 identifier = "%s"%(x.idNumber))\
                               for x in self.subsetStructures(requiredAttributes, sortAttribute=sortAttribute)],
                              fileName)

            
class Attribute():
    ''' Advantage of having a class over dict is that names of things standardized '''

    def __init__(self,
                 name     = None,
                 pointer  = None,
                 value    = None,
                 minValue = None,
                 maxValue = None,
                 notValue = None):

        # beware giving notValue=None - better to not give anything and have hasattr (do this another time)
        self.name     = name
        self.pointer  = pointer
        self.value    = value
        self.minValue = minValue
        self.maxValue = maxValue
        self.notValue = notValue
