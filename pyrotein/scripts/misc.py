#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from sklearn.neighbors import NearestNeighbors


def nearest_nbr(atm, nbrs):
    ''' Find the nearest neighbour for atm from nbrs.  

        Why not use `X = np.concatenate( (atm[np.newaxis,:], nbrs) )`???

        Assue user pass in Python list only.  DON'T recommend user to pass in
        numpy array because there is no benefit but hurdle that comes with
        numpy.  
    '''
    # Form a data structure for NN algorithm...
    X = np.concatenate( ([atm], nbrs) )

    # The NN procedure...
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(X)
    _, ind = nbrs.kneighbors(X)

    # Fetch the index of nearest neighbour...
    nbr_ind = ind[0][1]

    return X[nbr_ind]
