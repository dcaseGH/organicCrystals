import os
import shutil

---------------   Copy certain cryOpt files to PCM  --------------

eCut = -100.
#coDir  = '/home/dcase/FlufenamicAcid/CrystalPredictorFullRun/cryOpt/unique_pool'
pcmDir = '/home/dcase/FlufenamicAcid/CrystalPredictorFullRun/cryOpt/PCM/unique_pool'
runningList = []
print len(ffaCSPData.subsetStructures(requiredAttributes=[Attribute(pointer= 'crystalOptimiserComplete',
                                                                    value  = True),
                                                          Attribute(pointer={'crystalOptimiserData': 'finalStructureFileExists'},
                                                                    value=True)]))
for f in ffaCSPData.subsetStructures(requiredAttributes=[Attribute(pointer= 'crystalOptimiserComplete',
                                                                   value  = True),
                                                         Attribute(pointer={'crystalOptimiserData': 'finalStructureFileExists'},
                                                                   value=True)]):
    if f.getAttribute({'crystalOptimiserData': 'eTotalFinal'}) < eCut:
        if not os.path.isdir(pcmDir + '/' + str(f.idNumber) + '-1'):
            os.makedirs(pcmDir + '/' + str(f.idNumber) + '-1')
            shutil.copy(f.getAttribute({'crystalOptimiserData': 'finalStructurePointer'}),
                        pcmDir + '/' + str(f.idNumber) + '-1/starting.res')
            runningList.append(f.idNumber)

with open(pcmDir.replace('unique_pool', 'list'), 'w') as outf:
    outf.write('\n'.join(map(str, runningList+savedRunningList)))


-----------------------------------------------------------------------
