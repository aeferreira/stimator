"""S-timator : DEMO of dX/dt solution."""
from stimator import *

print __doc__
print
print """Vectors of dX/dt are computed by transformation:

transformation = m.getdXdt()
sdxdt.apply_transf(transformation)

-----------------------------------------------------------
"""

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

def step (t, t_stimulus, B):
    if t < t_stimulus:
        return 0.0
    else:
        return B
m.Bstep = forcing(step)

transformation = m.getdXdt()
svars = solve(m, tf = 6.0, npoints = 1000, title = 'X')
sdxdt = solve(m, tf = 6.0, npoints = 5000, title = 'dX/dt').apply_transf(transformation)

plot([svars,sdxdt], show = True)
