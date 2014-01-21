"""S-timator : DEMO of analysis module."""
from stimator import *
from stimator import __version__
from time import time
from numpy import append, linspace

print __doc__
print "S-timator version", __version__.fullversion, "(%s)" %__version__.date
print
print '---------------- EXAMPLE 1 ------------------'
m1 = read_model(examples.models.glyoxalases.text)
print (m1)

solution1 = solve(m1, tf=4030.0)

print '--- Last time point ----'
print 'At t =', solution1.t[-1]
for x, value in solution1.last:
    print "%-8s= %f" % (x, value)

#print '---------------- EXAMPLE 2 ------------------'
m2 = read_model(examples.models.branched.text)

#print m2
times = append(linspace(0.0, 5.0, 500), linspace(5.0, 10.0, 500))

solution2 = solve(m2, tf=10.0, times=times)

#print '---------------- EXAMPLE 3 ------------------'
m3 = read_model(examples.models.ca.text)
#print m3

solution3 = solve(m3, tf=8.0, npoints=2000)

print '---------------- EXAMPLE 4 ------------------'
m4 = read_model(examples.models.rossler.text)
print m4

solution4 = solve(m4, tf=100.0, npoints=2000, outputs="x1 x2 x3")


def transformation(vars, t):
    if t > 40.0:
        return (vars[0] - 5.0, vars[1], vars[2])
    else:
        return (-5.0, vars[1], vars[2])

solution4.apply_transf(transformation)

plot([solution1, solution2, solution3, solution4])

print '---------------- EXAMPLE 5 ------------------'
m5 = read_model(examples.models.lorentz.text)
print m5['title']
s = Solutions(m5['title'])
ms = ModelSolver(m5, tf=20.0,
                     npoints=20000,
                     outputs="x",
                     changing_pars="init.x")

for stimulus in 1.0, 1.01, 1.02:
    s += ms.solve(title='x(0) = %g' % stimulus, par_values=[stimulus])
plot(s, superimpose=True)

print '---------------- EXAMPLE 6 ------------------'
scantitle = "CICR model: Effect of stimulus on citosolic calcium"
m = read_model("""
title Calcium Spikes
v0         = -> Ca, 1
v1         = -> Ca, k1*B*step(t, 1.0)
k1         = 7.3
B          = 0.4
export     = Ca ->  , 10 ..
leak       = CaComp -> Ca, 1 ..

v2         = Ca -> CaComp, \
                  65 * Ca**2 / (1+Ca**2)
v3         = CaComp -> Ca, \
                  500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)
init       = state(Ca = 0.1, CaComp = 0.63655)
""")

s = Solutions(scantitle)
## print m


def mytransformation(B, Ca, CaComp):
    return B, Ca, CaComp
mytransformation.names = "stimulus", "Ca cit", 'Ca comp'

time0 = time()
print 'starting'
ms = ModelSolver(m, tf=6.0,
                    npoints=10000,
                    outputs="Ca CaComp",
                    changing_pars="B")
for stimulus in 0.0, 0.2, 0.4, 0.78:
    s += ms.solve(title='stimulus = %g' % stimulus, par_values=[stimulus])

## for stimulus in 0.0, 0.2, 0.4, 0.78:
##     m.B = stimulus
##     s += solve(m, tf = 6.0, npoints = 10000,
##     title = 'stimulus = %g'% (m.B), outputs="Ca CaComp")#mytransformation)

print 'done in', time() - time0, 's'
plot(s, ynormalize=True, show=True)
#plot(s, superimpose=True)
