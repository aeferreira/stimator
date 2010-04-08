#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

from model import *
from analysis import *

#print '---------------- EXAMPLE: CICR model ------------------'
def step (t, t_stimulus, B):
    if t < t_stimulus:
        return 0.0
    else:
        return B

m = Model("Calcium Spikes")
m.v0 = react(" -> Ca", 1)
m.v1 = react(" -> Ca", rate = "Bstep*k1")
m.Bstep = forcing(step)
m.k1 = 7.3
m.B  = 0.4
m.t_stimulus = 1.0
m.export = react(" Ca -> ", 10)
m.leak   = react("CaComp -> Ca", 1)
m.v2     = react("Ca -> CaComp", rate = "65 * Ca**2 / (1+Ca**2)")
m.v3     = react("CaComp -> Ca", rate = "500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)")
m.init   = state(Ca = 0.1, CaComp = 0.63655)

m.B = 0.4
svars = solve(m, tf = 6.0, npoints = 1000, title = 'X')
sdxdt = solve(m, tf = 6.0, npoints = 5000, title = 'dX/dt') #TODO: solutions should be 'clonable'.

transformation = m.getdXdt()
sdxdt.apply_transf(transformation)

plot([svars,sdxdt])
