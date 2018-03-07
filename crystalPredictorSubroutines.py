def analyseEditedLog(logFile):
    ''' point at something called Analyse_edited.log after the standard postprocessing scripts
        returns a list of dictionaries with information about the structures (energy ranked) '''
    with open(logFile, 'r') as readF:

        return[{'clusteredStructures': int(x.split('structures')[0].split('has')[1]),
                'Utot'               : float(x.split('Utot:')[1].split('Density')[0]),
                'density'            : float(x.split("\n")[2].split('Density:')[1]),
                'spaceGroup'         : x.split("\n")[3].split()[0]}
                for x in readF.read().split('Cluster:')[1:]]
