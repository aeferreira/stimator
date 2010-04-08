#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
import sys
import os.path

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir))

from stimator import *
from stimator.analysis import *
from time import time

#print '---------------- EXAMPLE: CICR model ------------------'
#~ def step (t):
    #~ if t < m.t_stimulus:
        #~ return 0.0
    #~ else:
        #~ return m.B

def step (t,B, t_stimulus):
    if t < t_stimulus:
        return 0.0
    else:
        return B

m = Model("Calcium Spikes")
m.v0 = " -> Ca", 1
m.v1 = " -> Ca",  "Bstep*k1"
m.Bstep = forcing(step)
m.k1 = 7.3
m.B  = 0.4
m.t_stimulus = 1.0
m.export = "Ca -> ", 10
m.leak   = "CaComp -> Ca", 1
m.v2     = "Ca -> CaComp",          "65 * Ca**2 / (1+Ca**2)"
m.v3     = ("CaComp -> Ca", 
            "500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)")
m.init   = state(Ca = 0.1, CaComp = 0.63655)

modeldef="""
title: Calcium Spikes
v0 : -> Ca, 1
v1 : -> Ca, rate = Bstep*k1
Bstep = forcing(step)
k1 = 7.3
B  = 0.4
t_stimulus = 1.0
export : Ca ->       , 10
leak   : CaComp -> Ca, 1
v2     : Ca -> CaComp, 65 * Ca**2 / (1+Ca**2)
v3     : CaComp -> Ca, 500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)
init   : state(Ca = 0.1, CaComp = 0.63655)

"""

#print m
s = Solutions("CICR model: Effect of stimulus on citosolic calcium")

def mytransformation(B, Ca, CaComp):
    return B, Ca, CaComp
mytransformation.names = "stimulus", "Ca cit", 'Ca comp'

time0 = time()
print 'starting'
for stimulus in 0.0, 0.2, 0.4, 0.78:
    m.B = stimulus
    s += solve(m, tf = 6.0, npoints = 1000, title = 'stimulus = %g'% (m.B), outputs="Ca CaComp")#mytransformation)

print 'done in', time()-time0, 's'
plot(s,ynormalize = True)
#plot(s, superimpose=True)
