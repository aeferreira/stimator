#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""
import sys
import os.path

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from numpy import *
from stimator import *
from stimator.fim import *


def test():
    
    print '---------------- EXAMPLE 1 ------------------'
    m = Model("Glyoxalase system")
    m.glo1 = react("HTA -> SDLTSH", rate = "V1*HTA/(Km1 + HTA)")
    m.glo2 = react("SDLTSH -> "   , rate = "V2*SDLTSH/(Km2 + SDLTSH)")
    m.V1  = 2.57594e-05
    m.Km1 = 0.252531
    m.V2  = 2.23416e-05
    m.Km2 = 0.0980973
    m.init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)

    #~ print m

    #~ solution1 = solve(m1, tf = 4030.0)

    #~ #print '---------------- EXAMPLE 3 ------------------'
    #~ m = Model("Calcium Spikes")
    #~ m.v0 = " -> Ca", 1
    #~ m.v1 = react(" -> Ca", rate = "B*k1")
    #~ m.k1 = 7.3
    #~ m.B  = 0.4
    #~ m.export = " Ca -> ", 10
    #~ m.leak   = "CaComp -> Ca", "kleak * CaComp"
    #~ m.kleak  = 1.0
    #~ m.v2     = "Ca -> CaComp", "65 * Ca**2 / (1+Ca**2)"
    #~ m.v3     = "CaComp -> Ca", "V3*CaComp**2/(CaComp**2+K3) * Ca**4/(Ca**4 + 0.6561)"
    #~ m.K3     = 4
    #~ m.V3     = 500
    #~ m.init   = state(Ca = 0.1, CaComp = 0.63655)

    #print m

    
    print '\nAdding sensitivity ODEs -------------------------'
    #pars = ["B", "k1", "K3"] # for calcium model
    pars = "V1 Km1".split()
    npars = len(pars)
    print 'npars =', npars
    nvars = len(m.variables)
    print 'nvars =', nvars
    nsens = npars * nvars
    print 'nsens =', nsens
    
    add_dSdt_to_model(m, pars)
    #print m
    
    #THIS WILL BE PART OF FIM    
    print '\nSensitivities at t = 3000 ----------------------'
    sol = solve(m, tf = 4030.0)

    i3000 = sol.i_time(3000)
    print 'i = ',i3000

    h = (sol.t[1] - sol.t[0])
    
    svec = sol.data[-nsens : , i3000]
    for i,x in enumerate(m.variables[-nsens:]):
        print x.name, svec[i]

    print '\nSensitivities at t = 3000 as a matrix ----------'
    S = matrix(svec.reshape((nvars, npars)))
    print S

    error_x = 0.001
    errorvar = error_x**2
    print '\nvar-covar matrix of measurement errors (const %g error) -------'% error_x
    #~ xvec = sol.data[ : nvars, i3000]
    #~ for i,x in enumerate(m.variables[: nvars]):
        #~ print x.name, xvec[i]
    MV = matrix(diag([errorvar,errorvar]))
    MVINV = linalg.inv(MV)
    print MVINV

    print '\nContribution to FIM at t = 3000 ----------------------'
    print 'h =', h
    ST = matrix(S.T)
    
    print "FIM (3000) ="
    FIMpoint = h * ST * MVINV * S
    print FIMpoint
    
    print '\nContribution to FIM at t = 3000 scaling by pars ------'
    parvalues = [getattr(m, p) for p in pars]
    #~ print 'parvalues =', parvalues
    P = matrix(diag(parvalues))
    print 'P ='
    print P
    S = S * P
    ST = matrix(S.T)
    print "FIM' (3000) ="
    FIMpoint = h * ST * MVINV * S
    print FIMpoint
    
    print '\nWhole timecourse ------------------------------------'
    ntimes = sol.data.shape[1]
    FIM = zeros((npars,npars))
    for i in range(ntimes):
        # S matrix
        svec = sol.data[-nsens : , i]
        S = matrix(svec.reshape((nvars, npars)))
        S = S * P

        #compute contribution at point i
        ST = matrix(S.T)
        FIMpoint = h * ST * MVINV * S
        
        #add to INVFIM
        FIM += FIMpoint
    print "FIM' ="
    print FIM
    PINV = linalg.inv(P)
    realFIM = PINV * FIM * PINV
    print "FIM ="
    print realFIM
    INVFIM = linalg.inv(FIM)
    INVFIM = P * INVFIM * P
    print "\nFIM-1 ="
    print INVFIM
    print "\nFIM-1 * FIM ="
    print dot(INVFIM, realFIM)
    print '\nParameters ------------------------------------'
    for i in range(len((pars))):
        print "%7s = %.5e +- %.5e" %(pars[i], parvalues[i], INVFIM[i,i]**0.5)
    
    print '--------------------------------------------'
    print '| CONSIDERING ONLY SDLTSH TIMECOURSE       |'
    print '--------------------------------------------'
    
    pars = "V1 Km1".split()
    npars = len(pars)
    vars = "SDLTSH".split()
    nvars = len(vars)
    nsens = npars * nvars
    
  
    print '\nSensitivities at t = 3000 ----------------------'
    indexes = []
    for x in vars:
        for i,y in enumerate(m.variables):
            dnames = y.name.split('_')
            if len(dnames) < 3: continue
            if dnames[1] == x and dnames[0] == 'd' and dnames[2] == 'd':
                indexes.append(i)

    indexes = array(indexes)
    svec = sol.data[indexes , i3000]
    for i,value in zip(indexes, svec):
        print m.variables[i].name, value

    print '\nSensitivities at t = 3000 as a matrix ----------'
    S = matrix(svec.reshape((nvars, npars)))
    print S

    error_x = 0.001
    errorvar = error_x**2
    print '\nvar-covar matrix of measurement errors (const %g error) -------'% error_x
    MV = matrix(diag([errorvar]))
    MVINV = linalg.inv(MV)
    print MVINV

    print '\nContribution to FIM at t = 3000 ----------------------'
    print 'h =', h
    ST = matrix(S.T)
    
    print "FIM (3000) ="
    FIMpoint = h * (ST * MVINV * S)
    print FIMpoint
    
    print '\nContribution to FIM at t = 3000 scaling by pars ------'
    parvalues = [getattr(m, p) for p in pars]
    # print 'parvalues =', parvalues
    P = matrix(diag(parvalues))
    print 'P ='
    print P
    S = S * P
    ST = matrix(S.T)
    print "FIM' (3000) ="
    FIMpoint = h * ST * MVINV * S
    print FIMpoint
    
    print '\nWhole timecourse ------------------------------------'
    ntimes = sol.data.shape[1]
    FIM = zeros((npars,npars))
    for i in range(ntimes):
        # S matrix
        svec = sol.data[indexes , i]
        S = matrix(svec.reshape((nvars, npars)))
        S = S * P

        #compute contribution at point i
        ST = matrix(S.T)
        FIMpoint = h * ST * MVINV * S
        
        #add to INVFIM
        FIM += FIMpoint
    print "FIM' ="
    print FIM
    PINV = linalg.inv(P)
    realFIM = PINV * FIM * PINV
    print "FIM ="
    print realFIM
    INVFIM = linalg.inv(FIM)
    INVFIM = P * INVFIM * P
    print "\nFIM-1 ="
    print INVFIM
    print "\nFIM-1 * FIM ="
    print dot(INVFIM, realFIM)
    print '\nParameters ------------------------------------'
    for i in range(len((pars))):
        print "%7s = %.5e +- %.5e" %(pars[i], parvalues[i], INVFIM[i,i]**0.5)
    
    print
    print '--------------------------------------------'
    print '|Checking results of computeFIM function   |'
    print '--------------------------------------------'
    #~ pars = "V1 Km1".split()
    #~ parvalues = [getattr(m, p) for p in pars]
    parsdict = dict (zip(pars, parvalues))

    print '\n--------- Example 1 ------------------'
    FIM1, invFIM1 = computeFIM4SingleTC(m, parsdict, "HTA SDLTSH".split(), 4030.0, 0.001)
    print "FIM ="
    print FIM1
    print 'FIM-1 ='
    print invFIM1
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    print '\n--------- Example 2 ------------------'
    FIM1, invFIM1 = computeFIM4SingleTC(m, parsdict, "SDLTSH".split(), 4030.0, 0.001)
    print "FIM ="
    print FIM1
    print 'FIM-1 ='
    print invFIM1
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    


if __name__ == "__main__":
    test()

