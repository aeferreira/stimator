#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
from stimator import *
from stimator.fim import add_dSdt_to_model


def test():
    
    print 'EXAMPLE 1 ------------------'
    print
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

    #~ #print ' EXAMPLE 2 ------------------'
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
    
    plot(plots)


if __name__ == "__main__":
    test()

