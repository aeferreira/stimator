"""S-timator : DEMO of analysis module."""
from stimator import *
from stimator import __version__
from time import time, sleep
from numpy import append, linspace

print __doc__
print "S-timator version", __version__.fullversion, "(%s)"%__version__.date
print
print '---------------- EXAMPLE 1 ------------------'
m1 = read_model("""
title Glyoxalase system
glo1 = HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
glo2 = SDLTSH ->    , V2*SDLTSH/(Km2 + SDLTSH)
V1  = 2.57594e-05
Km1 = 0.252531
V2  = 2.23416e-05
Km2 = 0.0980973
init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)
""")

print m1

solution1 = solve(m1, tf = 4030.0)

print '--- Last time point ----'
print 'At t =', solution1.t[-1]
for x,value in solution1.last:
    print "%-8s= %f" % (x, value)

#print '---------------- EXAMPLE 2 ------------------'
m2 = read_model("""
title Branched pathway
v1 = A -> B, k1*A
k1 = 10
v2 = B -> C, k2*B**0.5
k2 = 5
v3 = C -> D, k3*C**0.5
k3 = 2
v4 = C -> E, k4*C**0.5
k4 = 8
v5 = D ->  , k5*D**0.5
k5 = 1.25
v6 = E ->  , k6*E**0.5
k6 = 5
A  = 0.5
init = state(B = 2, C = 0.25, D = 0.64, E = 0.64)
""")

#print m2
times = append(linspace(0.0,5.0,500),linspace(5.0,10.0, 500))

solution2 = solve(m2, tf = 10.0, times=times)

#print '---------------- EXAMPLE 3 ------------------'
m3 = read_model("""
title Calcium Spikes
v0 =  -> Ca, 1
v1 =  -> Ca, B*k1
k1 = 7.3
B  = 0.4
export = Ca -> ,       10 ..
leak   = CaComp -> Ca, 1 ..
v2     = Ca -> CaComp, 65 * Ca**2 / (1+Ca**2)
v3     = CaComp -> Ca, 500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)
init   = state(Ca = 0.1, CaComp = 0.63655)
""")

#print m3

solution3 = solve(m3, tf = 8.0, npoints = 2000)

print '---------------- EXAMPLE 4 ------------------'
m4 = read_model("""
title Rossler
X1' = X2 - X3
X2' = 0.36 * X2 - X1
X3' = X1 *X3 - 22.5 * X3 - 49.6 * X1 + 1117.8
init = state(X1 = 19.0, X2 = 47, X3 = 50)
~ x3 = X3 -50.0
~ x1 = X1 -18.0
~ x2 = X2 -50.0
""")

print m4

solution4 = solve(m4, tf = 100.0, npoints = 2000, outputs="x1 x2 x3")

def transformation(vars,t):
    if t > 40.0:
        return (vars[0]-5.0, vars[1], vars[2])
    else:
        return (-5.0, vars[1], vars[2])

solution4.apply_transf(transformation)

sleep(1.0)

plot ([solution1, solution2, solution3, solution4])

s = Solutions("Rossler model: sensitivity to initial conditions")
ms = ModelSolver(m4,tf = 100.0, npoints = 2000, outputs="x1", changing_pars = "init.X1") 
for stimulus in 19.0, 19.02, 19.04:
    s += ms.solve(title = 'init.X1 = %g'% stimulus, par_values = [stimulus])
plot(s,superimpose=True)

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

s = Solutions("CICR model: Effect of stimulus on citosolic calcium")
print m
def mytransformation(B, Ca, CaComp):
    return B, Ca, CaComp
mytransformation.names = "stimulus", "Ca cit", 'Ca comp'

time0 = time()
print 'starting'
ms = ModelSolver(m,tf = 6.0, npoints = 10000, outputs="Ca CaComp", changing_pars = "B") 
for stimulus in 0.0, 0.2, 0.4, 0.78:
    s += ms.solve(title = 'stimulus = %g'% stimulus, par_values = [stimulus])

## for stimulus in 0.0, 0.2, 0.4, 0.78:
##     m.B = stimulus
##     s += solve(m, tf = 6.0, npoints = 10000, title = 'stimulus = %g'% (m.B), outputs="Ca CaComp")#mytransformation)

print 'done in', time()-time0, 's'
plot(s,ynormalize = True, show = True)
#plot(s, superimpose=True)
