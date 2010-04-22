#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *

def getCriteriumFunction(weights, tc):
    """Returns a function to compute the objective function (for each timecourse).
    
    the function has signature
    criterium(Y,i)
    Y is the predicted timecourse, for a given set of parameters.
    i is the index of the timecourse.
    The function returns a float.
    
    ydata is a list of ('experimental') timecourse data, 
    an array with shape (nt, nvars).
    
    weights can be:
    
    None         : no weighting (simple least squares, S = sum((Ypred-Yexp)**2))
    all others are weighted least squares, S = (Ypred-Yexp).T * W * (Ypred-Yexp)
    'demo'       : demo weighting  W = 1/j with j = 1,...,nvars
    
    
    """
    
    #mask series with NaN values.
    allvarindexes = []
    for data in tc:
        nt = data.ntimes
        varindexes = []

        for ivar in range(len(data.data)):
            #count NaN
            yexp = data[ivar]
            nnan = len(yexp[isnan(yexp)])
            if nnan >= nt-1: continue
            varindexes.append(ivar)
        allvarindexes.append(array(varindexes, int))
    

    if weights is None:
        def criterium(Y,i):
            d = (Y.T[allvarindexes[i]]- tc[i].data[allvarindexes[i]])
            return sum(d*d)
        return criterium

    if weights  == 'demo':
        W = []
        for i in range(len(tc)):
            W.append(array([1.0/(1+j) for j in range(allvarindexes[i])]))
        #print W
        def criterium(Y,i):
            d = (Y.T[allvarindexes[i]]- tc.data[i][allvarindexes[i]])
            return sum(d*W[i]*d)
        return criterium
    return None

def WSSD(differences, W):
    return sum(differences*W*differences)
