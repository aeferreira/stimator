"""S-timator : DEMO of analysis module."""
from stimator import *
from time import time

print __doc__

m = read_model("""
title Calcium Spikes
v0         = -> Ca, 1
v1         = -> Ca, Bstep*k1
Bstep      = 0.4
k1         = 7.3
B          = 0.4
t_stimulus = 1.0
export     = Ca ->  , 10 ..
leak       = CaComp -> Ca, 1 ..
    
v2         = Ca -> CaComp, \
                  65 * Ca**2 / (1+Ca**2)    
v3         = CaComp -> Ca, \
                  500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)
init       = state(Ca = 0.1, CaComp = 0.63655)
""")

def step (t,B, t_stimulus):
    if t < t_stimulus:
        return 0.0
    else:
        return B

m.Bstep = forcing(step)

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
