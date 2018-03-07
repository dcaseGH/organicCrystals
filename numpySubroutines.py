import numpy as np

def bestElementTwoMatricesAxiswise(m1, m2, bestM1 = 'max', bestM2 = 'max', axis1 = 0):
    ''' Returns a single element which is the best (max or min) in m1 along axis1
        in case of draw, returns best in m2
        in case of draw, returns lowest index in m2 '''

    if axis1 != 0:
        raiseException('implement other axis')

    if bestM1 == 'max':
        m1BestHits = np.max(m1, axis=axis1)
    elif bestM1 == 'min':
        m1BestHits = np.min(m1, axis=axis1)
    else:
        raiseException('bestM1 must be min or max')

    for i1 in xrange(m1BestHits.shape[0]):
        if bestM2 == 'min':
            m2BestHits = np.min([m2[i2, i1] for i2 in np.where(m1[:,i1] == m1BestHits[i1])[0]])
        elif bestM2 == 'max':
            m2BestHits = np.max([m2[i2, i1] for i2 in np.where(m1[:,i1] == m1BestHits[i1])[0]])
        else:
            raise Exception('bestM2 must be min or max')

        yield (np.where(m2[:, i1] == m2BestHits)[0][0], i1)
        

#use this as a testbe
#import numpy as np
#mat1 = np.array([[1,2,3],
#                 [31,2.3,3.2],
#                 [1,22,0.],
#                 [1.2,22,3444]])


#mat2 = np.array([[1,2,3],
#                 [31,2.3,3.2],
#                 [1,22,0.],
#                 [1.2, 0, 3444]])

#print [x for x in bestElementTwoMatricesAxiswise(mat1, mat2, bestM2='min')]
