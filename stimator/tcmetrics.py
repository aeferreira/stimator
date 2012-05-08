#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2010 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
from model import varnames


def _transform2array(vect):
    if isinstance(vect, float) or isinstance(vect, int):
        res = array((vect), dtype=float)
    elif isinstance(vect, list) or isinstance(vect, tuple):
        res = diag(array(vect, dtype=float))
    else:
        res = vect # is already an array (must be 2D)
    return res
    

def constError_func(vect):
    res = _transform2array(vect)
    def CE(x):
        return res
    return CE

def propError_func(vect):
    res = _transform2array(vect)
    def CE(x):
        return res * x
    return CE

def getFullTCvarIndexes(model, tcs):
    #mask series with NaN values.
    allmodelvarindexes, alltcvarindexes = [],[]
##     allvarindexes = []
    for data in tcs:
        nt = data.ntimes
        varindexes = []
        modelvarindexes = []

        for ivar in range(len(data.data)):
            #count NaN
            yexp = data[ivar]
            nnan = len(yexp[isnan(yexp)])
            if nnan >= nt-1: continue
            varindexes.append(ivar)
            vname = data.names[ivar]
            indx = varnames(model).index(vname)
            modelvarindexes.append(indx)
        alltcvarindexes.append(array(varindexes, int))
        allmodelvarindexes.append(array(modelvarindexes,int))
    print allmodelvarindexes, alltcvarindexes
    return allmodelvarindexes, alltcvarindexes

def getCommonFullVars(tcs):
    """Returns a list of names of variables that have full data in all timecourses."""
    common_names = []
    for itc,tc in enumerate(tcs):
        nt = tc.ntimes
        tcnames = tc.names
        for i,line in enumerate(tc.data):
            #count NaN
            yexp = line
            xname = tcnames[i]
            nnan = len(yexp[isnan(yexp)])
            if nnan >= nt-1:
                if xname in common_names:
                    common_names.remove(xname)
            else:
                if itc == 0:
                    common_names.append(xname)
    return common_names

def getRangeVars(tcs, varnames):
    ranges = [0.0 for i in range(len(varnames))]
    for ix,x in enumerate(varnames):
        for tc in tcs:
            yexp = tc[x]
            tpe = (max(yexp) - min(yexp))
            ranges[ix] = max(ranges[ix], tpe)
##             if tpe > ranges[ix]:
##                 ranges[ix] = tpe
    return ranges
    

def getCriteriumFunction(weights, model, tc):
    """Returns a function to compute the objective function (for each timecourse).
    
    the function has signature
    criterium(Y,i)
    Y is the predicted timecourse, for a given set of parameters.
    i is the index of the timecourse.
    The function returns a float.
    
    tc is a Solutions object holding ('experimental') timecourse data, 
    each timecourse has shape (nvars, ntimes).
    
    weights can be:
    
    None         : no weighting (simple least squares, S = sum((Ypred-Yexp)**2))
    all others are weighted least squares, S = (Ypred-Yexp).T * W * (Ypred-Yexp)
    'demo'       : demo weighting  W = 1/j with j = 1,...,nvars
    """
    
    allmodelvarindexes, alltcvarindexes = getFullTCvarIndexes(model,tc)    

    if weights is None:
        def criterium(Y,i):
            d = (Y.T[allmodelvarindexes[i]]- tc[i].data[alltcvarindexes[i]])
            return sum(d*d)
##             return sum(abs(d))
        return criterium

    if weights  == 'demo':
        W = []
        for i in range(len(tc)):
            W.append(array([1.0/(1+j) for j in range(alltcvarindexes[i])]))
        #print W
        def criterium(Y,i):
            d = (Y.T[allmodelvarindexes[i]]- tc.data[i][alltcvarindexes[i]])
            return sum(d*W[i]*d)
        return criterium
        
    ###TODO: weights not implemented
    return None
