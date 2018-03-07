import unittest

class TestInterfaces(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_ccdcToASE(self):
        from ioAndInterfaces import ccdcCrystalToASE
        from ccdc.io import CrystalReader
        ccdcCrystal = CrystalReader('FPAMCA16.res')[0]
        #print ccdcCrystal
        aseCrystal = ccdcCrystalToASE(ccdcCrystal)
        #print aseCrystal#.symbols
        # this is a Z'=4 30 atom molecule
        self.assertEqual(len(aseCrystal),
                         4 * 30)
        #test spacegroup and cell??

    def test_writingCIFFiles(self):
        pass
        
if __name__ == '__main__':
    unittest.main()
