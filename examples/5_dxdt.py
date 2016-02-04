"""S-timator : DEMO of dX/dt solution."""
from __future__ import print_function
from stimator import read_model, Solutions
from stimator.dynamics import getdXdt

print(__doc__, '\n')
print ("""Vectors of dX/dt are computed by transformation:

transformation = m.getdXdt()
sdxdt.apply_transf(transformation)

-----------------------------------------------------------
""")

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
init       : Ca = 0.1, CaComp = 0.63655
""")

x = m.solve(tf = 6.0, npoints = 5000, title = '$X$')
dxdt = x.copy(newtitle='$dX / dt$').apply_transf(getdXdt(m))

Solutions([x,dxdt]).plot(show = True)
