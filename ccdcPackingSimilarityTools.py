def returnNmatched_molecules(comparisonObject, crystal1, crystal2, returnRMSD=False):

    ''' All inputs are CCDC objects - this really just catches fails '''

    compareResult = comparisonObject.compare(crystal1, crystal2)

    if returnRMSD:
        try:
            return compareResult.nmatched_molecules, compareResult.rmsd
        except:
            return None, None
    else:
        try:
            return compareResult.nmatched_molecules
        except:
            return None
                
def matchStructureLists(list1,
                        list2,
                        allowMolecularDifferences=False,
                        returnRMSD   = False,
                        nRMSDMatches = 15):

    ''' Get number molecular matches for two lists of ccdc crystals 
        Return a numpy matrix '''

    from ccdc.crystal import PackingSimilarity
    import numpy as np
    
    ps = PackingSimilarity()
    ps.settings.allow_molecular_differences = allowMolecularDifferences
    ps.settings.packing_shell_size          = nRMSDMatches

    data = np.array([[returnNmatched_molecules(ps, cI, cJ, returnRMSD=returnRMSD)
                      for cI in list1]
                     for cJ in list2])
    if returnRMSD:
        return data[:,:,0], data[:,:,1]
    else:
        return data

#class MatchCSPStructureLists
