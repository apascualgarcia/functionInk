#!/usr/bin/env python

import pandas as pd
import numpy as np
import networkx as nx
import scipy
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import preprocessing
plt.ion()

# some simple spectral clustering test code
# def find_partition(evecs, k):
    # V = evecs[:,:k]
    # clust = KMeans(n_clusters = k)
    # clust.fit(V)
    # partition_vec = clust.labels_

    # return partition_vec

dataframe = pd.read_csv("./Rmrz89_Exhaust_Nest0.35_Conn0.276_CompConn0.15.long.txt",
                        sep="\t")

speciesA = dataframe.iloc[:,0].values
speciesB = dataframe.iloc[:,1].values
interactions = dataframe.iloc[:,2].values
num_edges = speciesA.shape[0]

G = nx.DiGraph()
for i in range(num_edges):
    G.add_edge(speciesA[i],speciesB[i],weight=interactions[i])

A = nx.adj_matrix(G)
Dabs = np.diag(abs(A).sum(axis=1).A.flatten())
Ls = Dabs -A

plt.figure()
plt.spy(A,markersize=1)


# some test code..
# variant 1 -- no diagonal scaling
# N1 = A*A.T
# N2 = A.T*A

# vals, vecs = scipy.sparse.linalg.eigsh(N1)
# vals2, vecs2 = scipy.sparse.linalg.eigsh(N2)


# variant 2 -- diagonal scaling
# S1 = Ls*Ls.T
# S2 = Ls.T*Ls
# vals3, vecs3 = scipy.sparse.linalg.eigsh(S1)
# vals4, vecs4 = scipy.sparse.linalg.eigsh(S2)

# u,s,v = scipy.sparse.linalg.svds(A,k=10)

# p0 = find_partition(vecs[::-1],2)
# p0 = find_partition(v.T[::-1],2)
# index  = np.argsort(p0)
# A2 = A[index,:]
# A2 = A2[:,index]
# plt.figure()
# plt.imshow(A2.A)
