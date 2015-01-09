"""S-timator : DEMO of  parameter scanning."""
from stimator import *
from stimator.dynamics import scan
from time import time, sleep
from numpy import linspace

def run_normal():
    print __doc__

    mdl = """
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
    init       : (Ca = 0.1, CaComp = 0.63655)
    """

    print mdl
    m = read_model(mdl)

    npoints = 1000

    title = "CICR model: Effect of stimulus on citosolic calcium"
    print (title)
    bvalues = {'B': (0.0, 0.1, 0.3, 0.5, 0.8, 1.0)}

    s = scan(m, bvalues, tf=10, npoints=npoints)
    s.plot(ynormalize = True, fig_size=(16,9), suptitlegend=title, show=True)
    #plot(s, superimpose=True)

if __name__ == "__main__":
    run_normal()

