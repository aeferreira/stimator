#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model related classes
# Copyright Ant�nio Ferreira 2006-2009
#----------------------------------------------------------------------------
from model import *
from numpy import *
from scipy import integrate

#import sys
#sys.path.append('..')
import pylab as p

def solve(model, tf = 1.0, npoints = 500, t0 = 0.0, initial = 'init'):

    salg=integrate._odepack.odeint
    m = model
    scale = 1.0
    f = m.scaled_dXdt(scale)
    y0 = copy(array(m.init))
    times = linspace (t0, tf, npoints)
    t  = copy(times)
    output = salg(f, y0, t, (), None, 0, -1, -1, 0, None, 
                    None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
    if output[-1] < 0: return (1.0E300)
    Y = output[0]
    #~ #append times to column 0
    #~ times = reshape(times, (-1,1))
    #~ return hstack((times,Y))
    return times, Y.T



print '---------------- EXAMPLE 1 ------------------'
m = Model("Glyoxalase system")
m.glo1 = react("HTA -> SDLTSH", rate = "V1*HTA/(Km1 + HTA)")
m.glo2 = react("SDLTSH -> "   , rate = "V2*SDLTSH/(Km2 + SDLTSH)")
m.V1  = 2.57594e-05
m.Km1 = 0.252531
m.V2  = 2.23416e-05
m.Km2 = 0.0980973
m.init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)

print m

t,solution = solve(m, tf = 4030.0)

#print t
#print solution
#print trf

print '--- Last time point ----'
print 'At t =', t[-1]
varnames = [x.name for x in m.variables]
for z in zip(varnames, solution[-1]):
    print "%-8s= %f" % z

#plot results...

f1 = p.figure(1)
p.subplot(221) 
for i, colour in enumerate(['r-', 'b-']):
    p.plot(t, solution[i], colour, label=m.variables[i].name)
p.grid()
p.legend(loc='best')
p.xlabel('')
p.ylabel('concentrations (mM)')
p.title(m.title)

print '---------------- EXAMPLE 2 ------------------'
m2 = Model("Branched pathway")
m2.v1 = react("A -> B", rate = "k1*A")
m2.k1 = 10
m2.v2 = react("B -> C", rate = "k2*B**0.5")
m2.k2 = 5
m2.v3 = react("C -> D", rate = "k3*C**0.5")
m2.k3 = 2
m2.v4 = react("C -> E", rate = "k4*C**0.5")
m2.k4 = 8
m2.v5 = react("D ->  ", rate = "k5*D**0.5")
m2.k5 = 1.25
m2.v6 = react("E ->  ", rate = "k6*E**0.5")
m2.k6 = 5
m2.A  = 0.5
m2.init = state(B = 2, C = 0.25, D = 0.64, E = 0.64)

print m2

t,solution = solve(m2, tf = 10.0)
#plot results...

p.subplot(222) 
for i, colour in enumerate(['r-', 'b-', 'g-', 'k-']):
    p.plot(t, solution[i], colour, label=m2.variables[i].name)
p.grid()
p.legend(loc='best')
p.xlabel('')
p.ylabel('X')
p.title(m2.title)

print '---------------- EXAMPLE 3 ------------------'
m3 = Model("Calcium Spikes")
m3.v0 = react(" -> Ca", rate = "k0")
m3.k0 = 1
m3.v1 = react(" -> Ca", rate = "B*k1")
m3.k1 = 7.3
m3.B  = 0.4
m3.export = react(" Ca -> ", rate = "k*Ca")
m3.k = 10
m3.leak = react("CaComp -> Ca", rate = "kf * CaComp")
m3.kf = 1
m3.v2 = react("Ca -> CaComp", rate = "65 * Ca**2 / (1+Ca**2)")
m3.v3 = react("CaComp -> Ca", rate = "500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)")
m3.init = state(Ca = 0.1, CaComp = 0.63655)

print m3

t,solution = solve(m3, tf = 8.0)
#plot results...

p.subplot(223) 
for i, colour in enumerate(['r-', 'b-']):
    p.plot(t, solution[i], colour, label=m3.variables[i].name)
p.grid()
p.legend(loc='best')
p.xlabel('')
p.ylabel('X')
p.title(m3.title)

print '---------------- EXAMPLE 4 ------------------'
m4 = Model("Rossler")
m4.v1 = react(" -> X1", rate = "X2 - X3")
m4.v2 = react(" -> X2", rate = "0.36 * X2 - X1")
m4.v3 = react(" -> X3", rate = "X1 *X3 - 22.5 * X3 - 49.6 * X1 + 1117.8")
m4.x3 = transf("X3 -50.0")
m4.x1 = transf("X1 -18.0")
m4.x2 = transf("X2 -50.0")
m4.init = state(X1 = 19.0, X2 = 47, X3 = 50)

print m4

t,solution = solve(m4, tf = 100.0, npoints = 1000)

#compute transf
f = m4.transf_func()
trf = apply_along_axis(f, 0, solution, 0.0)

#~ print len(t)
#~ print solution.shape
#~ print trf.shape

#plot results (transformations only)

p.subplot(224) 
for i, colour in enumerate(['r-', 'b-', 'g-']):
    p.plot(t, trf[i], colour, label=m4.transf[i].name)
p.grid()
p.legend(loc='best')
p.xlabel('')
p.ylabel('X')
p.title(m4.title)

p.show()