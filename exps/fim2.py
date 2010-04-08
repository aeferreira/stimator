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


def add_dSdt_to_model(m, pars):
    """Add sensitivity ODS to model, according to formula:
    
    dS/dt = df/dx * S + df/dp
    
    m is a model object
    pars are a list of parameter names
    """
    J = m.Jacobian_strings()
    dfdpstrs = m.dfdp_strings(pars)
    nvars = len(J)
    npars = len(pars)
    
    try:
        import sympy
    except:
        print 'ERROR: sympy module must be installed to generate sensitivity strings'
        raise
    _symbs = {}
    for x in m.variables:
        _symbs[x.name] = sympy.Symbol(x.name)
    for p in m.parameters:
        _symbs[p.name] = sympy.Symbol(p.name)
    for f in m.forcing:
        _symbs[f.name] = sympy.Symbol(f.name)
    
    #create symbols for sensitivities
    Smatrix = []
    for i in range(nvars):
        Smatrix.append([])
        for j in range(npars):
            Sname = "d_%s_d_%s" % (m.variables[i].name, pars[j])
            _symbs[Sname] = sympy.Symbol(Sname)
            Smatrix[i].append(Sname)

    for i in range(nvars):
        for j in range(npars):
            #compute string for dx/dt
            resstr = dfdpstrs[i][j]
            for k in range(nvars):
                resstr = resstr+ "+(%s)*(%s)"%(J[i][k], Smatrix[k][j])
            #~ print Smatrix[i][j]
            #~ print resstr
            #~ print type(resstr)
            _res = eval(resstr, None, _symbs)
            _dres = str(_res)
            if _dres == '0':
                _dres = '0.0'
            setattr(m, Smatrix[i][j], variable(_dres))
            setattr(m.init, Smatrix[i][j], 0.0)
    

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
    
    print '\nSolving with sensitivities ----------------------'
    plots = [solve(m, tf = 4030.0, title = 'X', outputs = "HTA SDLTSH")]
    for p in pars:
        plots.append(solve(m, tf = 4030.0, 
                              title = 'dX/d'+p , 
                              outputs = 'd_HTA_d_%s d_SDLTSH_d_%s'%(p,p)))
    #~ i1point2 = plots[0].i_time(1.2)
    #~ print i1point2
    
    #~ print '\nSolving with sensitivities ----------------------'
    #~ plots = [solve(m, tf = 3.0, npoints = 2700, title = 'X', outputs = "Ca CaComp")]
    #~ for p in pars:
        #~ plots.append(solve(m, tf = 3.0, 
                              #~ npoints = 2700, 
                              #~ title = 'dX/d'+p , 
                              #~ outputs = 'd_Ca_d_%s d_CaComp_d_%s'%(p,p)))
    #~ i1point2 = plots[0].i_time(1.2)
    #~ print i1point2
    
    #THIS WILL BE PART OF FIM    
    #~ print '\nSensitivities at t = 1.2 ----------------------'
    #~ sol = solve(m, tf = 3.0, npoints = 2700, title = 'all')
    #~ h = (sol.t[1] - sol.t[0])
    
    #~ svec = sol.data[-nsens : , i1point2]
    #~ for i,x in enumerate(m.variables[-nsens:]):
        #~ print x.name, svec[i]
    #~ print '\nSensitivities at t = 1.2 as a matrix ----------'
    #~ smatrix = matrix(svec.reshape((nvars, npars)))
    #~ print smatrix

    #~ print '\n diagonal 0.01 var-covar matrix of measurements ----------'
    #~ xvec = sol.data[ : nvars, i1point2]
    #~ for i,x in enumerate(m.variables[: nvars]):
        #~ print x.name, xvec[i]
    #~ MV = matrix(diag([0.01,0.01]))
    #~ MVINV = linalg.inv(MV)
    #~ print MVINV
    #~ print '\nContribution to FIM at t = 1.2 ----------------------'
    #~ ST = matrix(smatrix.T)
    
    #~ print "FIM (1.2) ="
    #~ FIM1point2 = h * ST * MVINV *smatrix
    #~ print FIM1point2
    
    #~ print '\nWhole timecourse ------------------------------------'
    #~ ntimes = sol.data.shape[1]
    #~ FIM = zeros((npars,npars))
    #~ for i in range(ntimes):
        #~ # S matrix
        #~ svec = sol.data[-nsens : , i]
        #~ smatrix = matrix(svec.reshape((nvars, npars)))

        #~ # 1% var-covar matrix of measurements
        #~ xvec = sol.data[ : nvars, i]
        #~ MV = matrix(diag([0.01,0.01]))
        
        #~ #conpute contribution at point i
        #~ ST = matrix(smatrix.T)
        #~ MVINV = linalg.inv(MV)
        #~ FIMpoint = h * ST * MVINV * smatrix
        
        #~ #add to INVFIM
        #~ FIM += FIMpoint
    #~ print "FIM ="
    #~ print FIM
    #~ INVFIM = linalg.inv(FIM)
    #~ print "\nFIM-1 ="
    #~ print INVFIM
    #~ print "\nFIM-1 * FIM ="
    #~ print dot(INVFIM, FIM)
    plot(plots)


if __name__ == "__main__":
    test()

