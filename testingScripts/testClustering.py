from cspTracking import SetCSPStructures, CSPStructure, generalCallErrorCatch, Attribute

tfaExpPCMData = SetCSPStructures.initFromCSV('/home/dcase/TolfenamicAcid/tolfenamicAcidDatabase.csv',
                                             calculations          ='pcm',
                                             selectIDNumberLambda  = lambda x: 'E' in x,
                                             idNumberApplyLambda   = lambda x: int(x.replace("E", "")),
                                             finalStructurePointer = lambda x: "/home/dcase/TolfenamicAcid/ExptData_PCM/KAXXAI_dfE%s_E%s.res"%(x,x),
                                             delim_whitespace      = False,
                                             useHeader             = True,
                                             specifiedStudyID      = 31,
                                             assertZPrime      = 1)
tfaExpPCMData.structures.append(CSPStructure(101,
                                             initPCMData = {'runDirectory': '/home/dcase/FlufenamicAcid/exptPolymorphs/TFA/PCM/unique_pool/4-1/',
                                                            'gaussianLog':  '/home/dcase/FlufenamicAcid/exptPolymorphs/TFA/PCM/unique_pool/5-1/mol1/act_mult1.log'}))
tfaExpPCMData.structures.append(CSPStructure(102,
                                             initPCMData = {'runDirectory': '/home/dcase/FlufenamicAcid/exptPolymorphs/TFA/PCM/unique_pool/5-1/',
                                                            'gaussianLog':  '/home/dcase/FlufenamicAcid/exptPolymorphs/TFA/PCM/unique_pool/5-1/mol1/act_mult1.log'}))

#tfaExpPCMData.makeCIF('expStructuresTFA.cif', sortAttribute = Attribute(pointer={'pcmData': 'eTotalFinal'}))
#62 matches 4 (after pcm minimisation - check this...)
assert( tfaExpPCMData.structureIDNumber(62).sameStructureAs(tfaExpPCMData.structureIDNumber(4), stage='pcmData') )
assert(len(tfaExpPCMData.structures) == 11)
tfaExpPCMData.clusterSet(stage = 'pcmData')#, clusterEnergyLimit = 10.)
print [(x.idNumber, x.pcmData['eTotalFinal'], generalCallErrorCatch(x, 'getAttribute', {'pcmData': 'clusterMatch'})) for x in tfaExpPCMData.structures]
assert(tfaExpPCMData.structureIDNumber(4).pcmData['clusterMatch'].idNumber == 62)

tfaExpPCMData.setGlobalAttributes([Attribute(pointer = {'pcmData': 'clusterMatch'},
                                             value   = None)])
assert(all([x.getAttribute({'pcmData': 'clusterMatch'}) == None for x in tfaExpPCMData.structures]))

tfaExpPCMData.setGlobalAttributes([Attribute(pointer = {'pcmData': 'clusterMatch'},
                                             value   = None)])
tfaExpPCMData.clusterSet(stage = 'pcmData', clusterEnergyLimit = 1.e-10)
assert(all([x.getAttribute({'pcmData': 'clusterMatch'}) == None for x in tfaExpPCMData.structures]))

tfaExpPCMData.setGlobalAttributes([Attribute(pointer = {'pcmData': 'clusterMatch'},
                                             value   = None)])
tfaExpPCMData.clusterSet(stage = 'pcmData', clusterEnergyTolerance = 1.e-10)
assert(all([x.getAttribute({'pcmData': 'clusterMatch'}) == None for x in tfaExpPCMData.structures]))
