import numpy as np
import unittest
import os

class TestTracking(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_csdTracker(self):
        import os
        from cspTracking import CSPStructure

        tempTracker = CSPStructure(1,
                                   initCrystalOptimiserData = 'exampleCOSummary.out')

        self.assertEqual(tempTracker.idNumber, 1)
        self.assertEqual(tempTracker.crystalOptimiserData['rhoFinal'],
                         1.4749)
        # assumes that Crys Opt was run where summary was found
        self.assertEqual(tempTracker.crystalOptimiserDir,
                         os.getcwd())

        #test getAttribute method
        self.assertEqual(tempTracker.getAttribute('idNumber'), 1)
        self.assertEqual(tempTracker.getAttribute({'crystalOptimiserData': 'rhoFinal'}),
                         1.4749)
        
        #test that it meets requested conditions (and "add" a space group for a check)
        tempTracker.crystalPredictorData = {'spaceGroup': None}
        from cspTracking import Attribute
        self.assertTrue(tempTracker.meetsConditions(requiredAttributes = [Attribute(pointer='idNumber'),
                                                                          Attribute(pointer={'crystalOptimiserData': 'rhoFinal'}),
                                                                          Attribute(pointer={'crystalPredictorData': 'spaceGroup'})]))

        #this example has a log file, but for a calculation which failed
        failTracker = CSPStructure(1,
                                   initCrystalOptimiserData = 'exampleCOFailDir/summary.out')
        self.assertFalse(failTracker.crystalOptimiserComplete)
        self.assertEqual(failTracker.crystalOptimiserError,
                         "electroStaticDBError_5833")

        #get energy per mol from cry opt and change it if there's a correction
        highZPStructure = CSPStructure(None,
                                       initCrystalOptimiserData = 'exampleHighZPCryOptSummary.out')
        self.assertEqual(highZPStructure.getAttribute({'crystalOptimiserData': 'zPrime'}),
                         4)
        np.testing.assert_approx_equal(highZPStructure.cryOptEnergyPerMol,
                                       -521.59936/4.)
        np.testing.assert_approx_equal(highZPStructure.getAttribute('cryOptEnergyPerMol'),
                                       -521.59936/4.)
        highZPStructure.cryOptEnergyPerMol += 1.
        np.testing.assert_approx_equal(highZPStructure.cryOptEnergyPerMol,
                                       1. + -521.59936/4.)

        # zero point energy from Jonas's AutoFree 1.0
        tempTracker.readFreeEnergyData('free_energy.log')
        self.assertEqual(tempTracker.getAttribute({'freeEnergyData': 'Evib'}),
                         -16.2758)
        self.assertEqual(tempTracker.getAttribute({'freeEnergyData': 'logPointer'}),
                         os.path.abspath('free_energy.log'))

        # pcm example with z'=3
        pcmData = CSPStructure(None, initPCMData = {'runDirectory': 'pcmExample',
                                                    'gaussianLog':  'pcmExample/mol*/*.log'})
        from cspTracking import HARTREE2KJPM
        self.assertAlmostEqual(pcmData.pcmData['eTotalFinal'],
                               (HARTREE2KJPM * (-1043.20905337 -1043.20886701 -1043.20928175) -428.1487)/3.)
        
    def test_SetCSPStructures(self):
        from cspTracking import SetCSPStructures

        tempSet = SetCSPStructures.initFromCrystalPredictorLog('exampleCP.log')

        self.assertEqual(len(tempSet.structures), 27414)
        self.assertEqual(tempSet.structures[0].crystalPredictorData['spaceGroup'], 'P-1')
        self.assertEqual(tempSet.structures[5].idNumber, 6)

        tempSet.structures[0].createCrystalOptimiserData('1-1/exampleCOSummary.out')
        self.assertEqual(tempSet.structures[0].crystalOptimiserData['rhoFinal'],
                         1.4749)

        #remove the data and add again with different method
        tempSet.structures[0].crystalOptimiserData = {}
        self.assertEqual(len(tempSet.structures[0].crystalOptimiserData.keys()), 0)
#        tempSet.addCrystalOptimiserInfoFromGlobString('1-1/*ummary.out')
        tempSet.addInfoFromGlobString('1-1/*ummary.out')
        self.assertEqual(tempSet.structures[0].crystalOptimiserData['rhoFinal'],
                         1.4749)
        self.assertTrue(tempSet.structures[0].crystalOptimiserComplete)

        # get free energy data from Jonas's log files
        tempSet.addFreeEnergyData(lambda x: 'free_energy.log')
        self.assertEqual(tempSet.structures[-1].getAttribute({'freeEnergyData': 'ZPE'}),
                         2.7239)

        
        #get a subset of data, and assert that it has a positive density
        from cspTracking import Attribute
        self.assertEqual(len(tempSet.subsetStructures([Attribute(pointer='crystalOptimiserComplete',
                                                                 value  = True)])),
                         len([x for x in tempSet.structures if x.crystalOptimiserComplete]))
        self.assertEqual(len(tempSet.subsetStructures([Attribute(pointer   = {'crystalOptimiserData': 'rhoFinal'},
                                                                 minValue  = 0.)])),
                         1)
        self.assertEqual(len(tempSet.subsetStructures([Attribute(pointer   = {'crystalOptimiserData': 'rhoFinal'},
                                                                 maxValue  = 0.)])),
                         0)
        self.assertEqual(len(tempSet.subsetStructures([Attribute(pointer   = {'crystalOptimiserData': 'rhoFinal'},
                                                                 minValue  = 0.1,
                                                                 maxValue  = 0.4)])),
                         0)
#        tempSet.addCrystalOptimiserInfoFromGlobString('2-1/*ummary.out')
        tempSet.addInfoFromGlobString('2-1/*ummary.out')
        self.assertEqual(len(tempSet.subsetStructures([Attribute(pointer   = {'crystalOptimiserData': 'rhoFinal'},
                                                                 minValue  = 1.1,
                                                                 maxValue  = 2.4)])),
                         2)
        #lowest one for rho has idNumber 2
        self.assertEqual([x.getAttribute('idNumber') for x in
                          tempSet.subsetStructures([Attribute(pointer   = {'crystalOptimiserData': 'rhoFinal'},
                                                              minValue  = 1.1,
                                                              maxValue  = 2.4)],
                                                   sortAttribute=Attribute(pointer   = {'crystalOptimiserData': 'rhoFinal'}))][0],
                         2)
                         
        #turn data into a pandas DataFrame
        dataTable = tempSet.toPandasDataFrame()
        self.assertEqual(dataTable['SG'][0], 'P-1')

        #test matching        
        from ccdc.io import CrystalReader
        from cspTracking import CSPStructure

        matchingSet = SetCSPStructures([CSPStructure(1, initCrystalOptimiserData='1-1/exampleCOSummary.out')])
        #match to itself just as an example
        self.assertTrue(matchingSet.structures[0].crystalOptimiserData['finalStructureFileExists'])
        matchDict = matchingSet.matchToStructureList([CrystalReader('1-1/lowest.res')[0]],
                                                     stage='crystalOptimiserData')
        np.testing.assert_array_equal(matchDict['matchesMatrix'],
                                      np.array([[15]]))
                                                 
        self.assertEqual(matchDict['candidateList'][0].idNumber,
                         1)

        #do same match but add RMSD data
        matchDict = matchingSet.matchToStructureList([CrystalReader('1-1/lowest.res')[0]],
                                                     returnRMSD=True,
                                                     stage='crystalOptimiserData')
        np.testing.assert_array_equal(matchDict['rmsdMatrix'],
                                      np.array([[0.]]))
        
        #test packing coefficient (what method does CCDC use??)
        self.assertAlmostEqual(matchingSet.structures[0].packingCoefficient(stage='crystalOptimiserData'),
                               0.69370122631)
        
        #get minimum energy of list
        self.assertEqual(matchingSet.minValueInSet(Attribute(pointer = {'crystalOptimiserData': 'eTotalFinal'})),
                         -133.44284)

        #test writing a CIF - tear down data at the end
        nameOfCIF = 'tempCIFWrite.cif'
        self.assertFalse(os.path.isfile(nameOfCIF))

        matchingSet.calculation = 'crystalOptimiser' # this should be known before matching really
        matchingSet.makeCIF(nameOfCIF)
        self.assertTrue(os.path.isfile(nameOfCIF))
        os.remove(nameOfCIF)        

        #initiate via CSV then match
        csvStructures = SetCSPStructures.initFromCSV('example.csv',
                                                     idNumberIndex     = 0,
                                                     finalEnergyIndex  = 1,
                                                     finalDensityIndex = 3,
                                                     finalStructureDir = os.path.dirname(os.path.realpath(__file__)),
                                                     assertZPrime      = 1.234)
        self.assertEqual(csvStructures.structures[0].idNumber,
                         978)
        self.assertEqual(csvStructures.structures[0].getAttribute({'crystalOptimiserData': 'zPrime'}),
                         1.234)
        self.assertAlmostEqual(csvStructures.structures[0].getAttribute({'crystalOptimiserData': 'eTotalFinal'}),
                               -139.35473)
        self.assertAlmostEqual(csvStructures.structures[0].getAttribute({'crystalOptimiserData': 'rhoFinal'}),
                               1.26410)
        print 'check match'
        
        #use general get and set attributes
        newStructureSet = SetCSPStructures([CSPStructure(i) for i in range(1,5)], calculation = 'crystalOptimiser')
        newStructureSet.structures[2].setAttribute(Attribute(pointer={'crystalOptimiserData': 'rhoFinal'},
                                                             value = 23.4))
        self.assertEqual(newStructureSet.structures[2].getAttribute({'crystalOptimiserData': 'rhoFinal'}),
                         23.4)
        newStructureSet.setGlobalAttributes([Attribute(pointer={'crystalOptimiserData': 'rhoFinal'},
                                                       value=1.0)])
        self.assertEqual(newStructureSet.structures[2].getAttribute({'crystalOptimiserData': 'rhoFinal'}),
                         1.)
        self.assertEqual(newStructureSet.calculation,
                         'crystalOptimiser')
#        print newStructureSet.structures[0].__dict__
        print newStructureSet.statusString()

        #read some PCM data
        newStructureSet.addInfoFromGlobString('1-1/PCM/1-1', calculation = 'pcm')
        self.assertAlmostEqual(newStructureSet.structures[0].pcmData['eInterFinal'],
                               -149.0522)
        self.assertAlmostEqual(newStructureSet.structures[0].pcmData['eIntraHartrees'],
                               -1043.20561970)
        self.assertTrue(newStructureSet.structures[0].pcmData['finalStructureFileExists'])
        
        def tearDown(self):
            print 'tearing down', os.path.isfile('tempCIFWrite.cif')
#            for deleteFiles in ['tempCIFWrite.cif']:
#                os.remove(deleteFiles)

if __name__ == '__main__':
    unittest.main()
