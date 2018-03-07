#this could be a useful way to split things up - but ccdc has functionality to get components anyway
# the sklearn thing- I don't know what it does. The scipy appears easy and useful

import numpy as np
import scipy
from sklearn.neighbors import kneighbors_graph

X = scipy.linalg.block_diag([[0, 1], [1, 0]], [[0, 1, 0],[1, 0, 1], [0, 1, 0]])
print X
#X[:,4], X[4, :], X[:, 1], X[1, :] = X[:,1], X[1, :], X[:, 4], X[4, :] 
#print X
A = kneighbors_graph(X, 2, mode='connectivity', include_self=True)

print A.toarray()
print A.__dict__

print 'sparse'
xSparse = scipy.sparse.csgraph.connected_components(X, connection='strong')
print xSparse

Y = np.array([[0, 0, 1, 0, 0],
              [0, 0, 0, 1, 0],
              [1, 0, 0, 0, 0],
              [0, 1, 0, 0, 1],
              [0, 0, 0, 1, 0]])

ySparse = scipy.sparse.csgraph.connected_components(Y, connection='strong')
print ySparse

Z = np.array([[0, 0, 1, 0, 0, 0, 0],
              [0, 0, 0, 1, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 1, 0, 0],
              [0, 0, 0, 1, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 1],
              [0, 0, 0, 0, 0, 1, 0]])

zSparse = scipy.sparse.csgraph.connected_components(Z, connection='strong')
print zSparse
