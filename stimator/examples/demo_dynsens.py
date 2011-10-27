#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : DEMO of dynamic sensitivities."""

from numpy import *
from stimator import *
from stimator.analysis import add_dSdt_to_model

print __doc__
print
print """Sensitivity ODEs are added to model, according to formula:

     dS/dt = df/dx * S + df/dp
-----------------------------------------------------------
"""


glos = """
title Glyoxalase system
glo1 : HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
glo2 : SDLTSH -> ,    V2*SDLTSH/(Km2 + SDLTSH)

V1  = 2.57594e-05
Km1 = 0.252531
V2  = 2.23416e-05
Km2 = 0.0980973

init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)
"""
m = read_model(glos)
print m

print '\nAdding sensitivity ODEs -------------------------'
#pars = ["B", "k1", "K3"] # for calcium model
pars = "V1 Km1".split()
npars = len(pars)
print 'npars =', npars
nvars = len(varnames(m))
print 'nvars =', nvars
nsens = npars * nvars
print 'nsens =', nsens

add_dSdt_to_model(m, pars)
#print m

print '\nSolving with sensitivities...'
sol = solve(m, tf = 4030.0)
plots = [sol.copy(names = "HTA SDLTSH", newtitle = 'X')]
for p in pars:
    plots.append(sol.copy(names = 'd_HTA_d_%s d_SDLTSH_d_%s'%(p,p), 
                          newtitle ='dX/d'+p))

print '\nDONE!'
plot(plots, show = True)
