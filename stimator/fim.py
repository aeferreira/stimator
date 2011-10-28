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
import expcov

sympy_installed = True
try:
    import sympy
except:
    print 'ERROR: sympy module must be installed to generate sensitivity strings'
    sympy_installed = False

def __computeNormalizedFIM(model, pars, vars, timecoursedata, expCOV):
    
    """Computes FIM normalized by parameter values.
    
    model is a model object.
    
    pars is a dicionary of parameter names as keys and parameter values as values,
    or list of (name,value) pairs
        names of the form init.<name> are possible.
        
    vars is a list of names of variables to be considered in the timecourses
    
    timecoursecollection is a Solutions object 
       (initial values are read from these)
    
    expCOV is a function wich computes the 
        experimental variance-covariance matrix of variables.
        function has signature f(x), where x is the numpy array of variables.
        NOTE: module expcov defines some 'stock' error functions.
    """
    m  = model.clone()
    
    #ensure m has init attr
    inits = {}
    for x in varnames(m):
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
##         print "****** m.init =",m.init
        sols += solve(m, tf = tc.t[-1])
    
    #compute P

    # P is the diagonal matrix of parameter values
    P = matrix(diag(array(parvalues, dtype=float)))
    
    #keep indexes of variables considered
    vnames = varnames(m)
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
##     print "INDEXES =", indexes
    indexes = array(indexes)

    tcFIM = []
    for sol in sols:
##         xvec = sol.data[indexes , -1]
##         print "######### LAST point"
##         print xvec
        #compute integral of ST * MVINV * S
        ntimes = len(sol.t)
        FIM = zeros((npars,npars))
        for i in range(1,ntimes):
            #time step
            h = (sol.t[i]-sol.t[i-1])
            
            #S matrix
            svec = sol.data[indexes , i]
            S = matrix(svec.reshape((nvars, npars)))
            S = S * P # scale with par values
            ST = matrix(S.T)

            # MVINV is the inverse of measurement COV matrix
            xvec = sol.data[xindexes , i]
            error_x = expCOV(xvec)
            MV = matrix(error_x**2)
            MVINV = linalg.inv(MV)
            
            #contribution at point i (rectangle's rule)
            FIMpoint = h * ST * MVINV * S
            FIM += FIMpoint
        tcFIM.append(FIM)
    
    sumFIM = sum(dstack(tcFIM), axis=2)
    
    return sumFIM, P

def computeFIM(model, pars, vars, TCs, COV):
    FIM, P = __computeNormalizedFIM(model, pars, vars, TCs, COV)
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
    errors = expcov.constError_func(errors)
    errorsSDLonly = 0.001
    errorsSDLonly = expcov.constError_func(errorsSDLonly)

    print '\n------------------------------------------------'
    print 'Glyoxalase model, 1 timecourse, parameters %s and %s'% (pars[0],pars[1])
    print 'Timecourse with HTA and SDLTSH'
    FIM1, invFIM1 = computeFIM(m, parsdict, "HTA SDLTSH".split(), sols, errors)
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    print '\n------------------------------------------------'
    print 'Glyoxalase model, 1 timecourse, parameters %s and %s'%(pars[0],pars[1])
    print 'Timecourse with SDLTSH only'
    FIM1, invFIM1 = computeFIM(m, parsdict, "SDLTSH".split(), sols, errorsSDLonly)
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
    
    FIM1, invFIM1 = computeFIM(m, parsdict, "HTA SDLTSH".split(), sols, errors)
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    
    print '\n------------------------------------------------'
    print 'Glyoxalase model, 2 timecourses, parameters %s and %s'% (pars[0],pars[1])
    print 'Timecourses with SDLTSH only'
    
    FIM1, invFIM1 = computeFIM(m, parsdict, ["SDLTSH"], sols, errorsSDLonly)
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)

if __name__ == "__main__":
    test()

