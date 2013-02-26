#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2010 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
from model import *
from dynamics import add_dSdt_to_model, dXdt_strings
from analysis import solve
import timecourse
import tcmetrics

sympy_installed = True
try:
    import sympy
except:
    print 'ERROR: sympy module must be installed to generate sensitivity strings'
    sympy_installed = False

def __computeNormalizedFIM(model, pars, timecoursedata, expCOV, vars = None):
    
    """Computes FIM normalized by parameter values.
    
    model is a model object.
    
    pars is a dicionary of parameter names as keys and parameter values as values,
    or list of (name,value) pairs
        names of the form init.<name> are possible.
        
    timecoursecollection is a Solutions object 
       (initial values are read from these)
    
    expCOV is a function wich computes the 
        experimental variance-covariance matrix of variables.
        function has signature f(x), where x is the numpy array of variables.
        NOTE: module expcov defines some 'stock' error functions.

    vars is a list of names of variables to be considered in the timecourses
    
    """
    m  = model.clone()
    
    #ensure m has init attr
    inits = {}
    for x in m().varnames:
        inits[str(x)] = 0.0
    setattr(m, 'init', state(**inits))
    
    if isinstance(pars,dict):
        parnames = pars.keys()
        parvalues = pars.values()
    else:
        parnames = [n for (n,v) in pars]
        parvalues = [v for (n,v) in pars]
    
    for n,v in zip(parnames,parvalues):
        obj = m
        names = n.split('.')
        name = names[-1]
        for n in names[:-1]:
            obj = getattr(obj, n)
        setattr(obj, name, v)
    
    # Adding sensitivity ODEs
    npars = len(pars)
    nvars = len(vars)
    Snames = add_dSdt_to_model(m, parnames)
##     print '\ndXdt_strings(): --------------------------'
##     for xname,dxdt in dXdt_strings(m):
##         print '(d %s /dt) ='%(xname),dxdt
##     print
##     print "*****Snames =", Snames

    sols = timecourse.Solutions()
    for tc in timecoursedata:
        #set init and solve
        for n,v in tc.state_at(0.0): #TODO: may not start at zero
            setattr(m.init, n, v)
        sols += solve(m, tf = tc.t[-1])
    
    #compute P

    # P is the diagonal matrix of parameter values
    P = matrix(diag(array(parvalues, dtype=float)))
    
    #keep indexes of variables considered
    if vars is not None:
        vnames = m().varnames
        xindexes = []
        for vname in vars:
            for i,y in enumerate(vnames):
                if y == vname:
                    xindexes.append(i)
        xindexes = array(xindexes)
        #compute indexes of sensitivities
        # search pattern d_<var name>_d_<parname> in variable names
        indexes = []
        for vname in vars:
            for Sname in Snames:
                if Sname[0] == vname:
                    indexes.append(vnames.index(Sname[2]))
        indexes = array(indexes)
        # insert indexes in sols as attributes
        for sol in sols:
            sol.xindexes = xindexes
            sol.indexes = indexes

    tcFIM = []
    for sol in sols:
        #compute integral of ST * MVINV * S
        ntimes = len(sol.t)
        FIM = zeros((npars,npars))
        for i in range(1,ntimes):
            #time step
            h = (sol.t[i]-sol.t[i-1])
            
            #S matrix
            svec = sol.data[sol.indexes , i]
            S = matrix(svec.reshape((nvars, npars)))
            S = S * P # scale with par values
            ST = matrix(S.T)

            # MVINV is the inverse of measurement COV matrix
            xvec = sol.data[sol.xindexes , i]
            error_x = expCOV(xvec)
            MV = matrix(error_x**2)
            MVINV = linalg.inv(MV)
            
            #contribution at point i (rectangle's rule)
            FIMpoint = h * ST * MVINV * S
            FIM += FIMpoint
        tcFIM.append(FIM)
    
    sumFIM = sum(dstack(tcFIM), axis=2)
    
    return sumFIM, P

def computeFIM(model, pars, TCs, COV, vars = None):
    FIM, P = __computeNormalizedFIM(model, pars, TCs, COV, vars)
    # compute inverse of P to "descale"
    PINV = linalg.inv(P)
    # compute FIM and its inverse (lower bounds for parameter COV matrix)
    realFIM = PINV * FIM * PINV
    INVFIM = linalg.inv(FIM)
    INVFIM = P * INVFIM * P
    check = dot(INVFIM, realFIM)
    #TODO: MAKE THIS CHECK USEFUL
    return realFIM, INVFIM

def test():
    
    m = Model("Glyoxalase system")
    m.glo1 = react("HTA -> SDLTSH", rate = "V*HTA/(Km1 + HTA)", pars=dict(V= 2.57594e-05))
    m.glo2 = react("SDLTSH -> "   , rate = "V2*SDLTSH/(Km2 + SDLTSH)")
    m.V   = 2.57594e-05
    m.Km1 = 0.252531
    m.V2  = 2.23416e-05
    m.Km2 = 0.0980973
    m.init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)
    
##     pars = "V Km1".split()
##     parvalues = [m.V, m.Km1]
    pars = "glo1.V Km1".split()
    parvalues = [m.glo1.V, m.Km1]
    
    sols = timecourse.Solutions()
    sols += solve(m, tf = 4030.0) 
    
    parsdict = dict (zip(pars, parvalues))
    
    errors = (0.01,0.001)
    errors = tcmetrics.constError_func(errors)
    errorsSDLonly = 0.001
    errorsSDLonly = tcmetrics.constError_func(errorsSDLonly)

    print '\n------------------------------------------------'
    print 'Glyoxalase model, 1 timecourse, parameters %s and %s'% (pars[0],pars[1])
    print 'Timecourse with HTA and SDLTSH'
    FIM1, invFIM1 = computeFIM(m, parsdict, sols, errors,"HTA SDLTSH".split())
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    print '\n------------------------------------------------'
    print 'Glyoxalase model, 1 timecourse, parameters %s and %s'%(pars[0],pars[1])
    print 'Timecourse with SDLTSH only'
    FIM1, invFIM1 = computeFIM(m, parsdict, sols, errorsSDLonly, "SDLTSH".split())
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    

    print '\n------------------------------------------------'
    print 'Glyoxalase model, 2 timecourses, parameters %s and %s'% (pars[0],pars[1])
    print 'Timecourses with HTA and SDLTSH'

    #generate 2nd timecourse
    m.init.SDLTSH = 0.001246154
    m.init.HTA = 0.2688
    sols += solve(m, tf = 5190.0)
    
    FIM1, invFIM1 = computeFIM(m, parsdict, sols, errors, "HTA SDLTSH".split())
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    
    print '\n------------------------------------------------'
    print 'Glyoxalase model, 2 timecourses, parameters %s and %s'% (pars[0],pars[1])
    print 'Timecourses with SDLTSH only'
    
    FIM1, invFIM1 = computeFIM(m, parsdict, sols, errorsSDLonly,["SDLTSH"])
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)

if __name__ == "__main__":
    test()

