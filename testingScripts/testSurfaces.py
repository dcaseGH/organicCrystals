import numpy as np
import unittest
import os

print '''Things to fix:
splitting molecules/using RDKit and losing decimal places
'''
class AnalyseClusters(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_makeSurface(self):
        from analyseClusters import getSpeciesListCIF, Surface, splitASEAtomsToMols

        structures = list(getSpeciesListCIF('CSPPCM_example.cif'))
        #take the z'=2 one - but hack it to make one atom per mol
        hackMols = splitASEAtomsToMols(structures[-2])
        from ase import Atoms

        mySurface = Surface(hackMols)

        #try simple cells and cell returning sensible cells
        _oldCell = mySurface.aseCell
        mySurface = Surface(hackMols, millerIndex = np.array([0,1,0]))
        np.testing.assert_array_almost_equal(np.dot(np.array([[0., 0., 1.],
                                                              [1., 0., 0.],
                                                              [0., 1., 0.]]),
                                                    _oldCell),
                                             mySurface.aseCell)
        mySurface = Surface(hackMols, millerIndex = np.array([1,0,0]))
        np.testing.assert_array_almost_equal(np.dot(np.array([[0., 0., 1.],
                                                              [1., 0., 0.],
                                                              [0., 1., 0.]]).T,
                                                    _oldCell),
                                             mySurface.aseCell)
        #this needs a cubic cell otherwise becomes big
        hackMols = [Atoms(symbols=['C'],  positions=[hackMols[0].get_center_of_mass()], cell=hackMols[0].cell, info=hackMols[0].info),
                    Atoms(symbols=['Xe'], positions=[hackMols[1].get_center_of_mass()], cell=hackMols[1].cell, info=hackMols[1].info)]
        hackMols = [Atoms(symbols=['Xe', 'C'], positions=[[0.25, 0.5, 0.75], [0.75, 1., 1.]], cell = [[1., 0., 0.], [ 0., 2., 0.],[0.,0.,3]], info=hackMols[0].info)]
        mySurface = Surface(hackMols, millerIndex = np.array([0,0,1]))
        np.testing.assert_array_almost_equal(np.diag([1., 2., 3.]),
                                             mySurface.aseCell)

        #try other cuts
        mySurface.writeASEAtoms('hack001.cif')
        mySurface = Surface(hackMols, millerIndex = np.array([0,1,1]))
        mySurface.writeASEAtoms('hack011.cif')
#        print mySurface.aseCell


    def test_printStuff(self):
        from analyseClusters import getSpeciesListCIF, Surface, splitASEAtomsToMols
        return
#        mySurface = Surface(splitASEAtomsToMols(structures[-2]), millerIndex = np.array([1,1,1]))
        structures = list(getSpeciesListCIF('CSPPCM_example.cif'))
        #take the z'=2 one - but hack it to make one atom per mol
        hackMols = splitASEAtomsToMols(structures[-2])
        from ase import Atoms
        #draw the cell (for some reason, it may help to see it)
        mySurface.asymmetricMolecules = mySurface.filledUnitCellMolecules()
        mySurface.writeASEAtoms('originalSetting.cif')
        _cifName = 'exampleSurface.cif'
        if os.path.isfile(_cifName):
            os.remove(_cifName)
        mySurface.writeASEAtoms(_cifName, pbc = False)
        self.assertTrue(os.path.isfile(_cifName))

    def test_CrystalSubroutines(self):
        from analyseClusters import getSpeciesListCIF, Crystal, splitASEAtomsToMols

        structures = list(getSpeciesListCIF('CSPPCM_example.cif'))
        #take the z'=2 one
        myCrystal = Crystal(splitASEAtomsToMols(structures[-2]))
        #spacegroup 7, z'2 = 4
        self.assertTrue(len(myCrystal.filledUnitCellMolecules()),
                        4)

        #check that in ase the convention is cartVec = fracVec.Cell 
        np.testing.assert_array_almost_equal(myCrystal.asymmetricMolecules[0].aseAtoms.get_center_of_mass(),
                                             np.dot(myCrystal.asymmetricMolecules[0].aseAtoms.get_center_of_mass(scaled=True),
                                                    myCrystal.aseCell))

        #check that in ase the convention is fracVec = cartVec.(Cell^-1) 
        np.testing.assert_array_almost_equal(myCrystal.asymmetricMolecules[0].aseAtoms.get_center_of_mass(scaled=True),
                                             np.dot(myCrystal.asymmetricMolecules[0].aseAtoms.get_center_of_mass(),
                                                    np.linalg.inv(myCrystal.aseCell)))

    def test_surfaceDefinedLattice(self):
        print 'this really needs to be tested'


if __name__ == '__main__':
    unittest.main()
