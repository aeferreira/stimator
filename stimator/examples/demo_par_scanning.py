"""S-timator : DEMO of scanning parameter functions."""
from stimator import *
from time import time, sleep
from numpy import append, linspace

print __doc__

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

print m

def run_normal():
    s = Solutions(title="CICR model: Effect of stimulus on citosolic calcium")
    time0 = time()
    ms = ModelSolver(m,tf = 6.0, npoints = 1000, outputs="Ca CaComp", changing_pars = "B") 
    print 'starting'
    for stimulus in linspace(0.0,1.0,200):
        s += ms.solve(title = 'stimulus = %g'% stimulus, par_values = [stimulus])

    print 'using ModelSolver done in', time()-time0, 's'

    s = Solutions(title="CICR model: Effect of stimulus on citosolic calcium")

    time0 = time()
    for stimulus in linspace(0.0,1.0,200):
        m.B = stimulus
        s += solve(m, tf = 6.0, npoints = 1000, title = 'stimulus = %g'% (m.B), outputs="Ca CaComp")#mytransformation)

    print 'using solve done in', time()-time0, 's'

    s = Solutions(title="CICR model: Effect of stimulus on citosolic calcium")

    ms = ModelSolver(m,tf = 6.0, npoints = 1000, outputs="Ca CaComp", changing_pars = "B") 
    for stimulus in 0.0, 0.2, 0.4, 0.78:
        s += ms.solve(title = 'stimulus = %g'% stimulus, par_values = [stimulus])

    for stimulus in 0.0, 0.2, 0.4, 0.78:
        m.B = stimulus
        s += solve(m, tf = 6.0, npoints = 1000, title = 'stimulus = %g'% (m.B), outputs="Ca CaComp")#mytransformation)

    s.plot(ynormalize = True, show = True)
    #plot(s, superimpose=True)

def test():
    s = Solutions("CICR model: Effect of stimulus on citosolic calcium")
    ms = ModelSolver(m,tf = 6.0, npoints = 1000, outputs="Ca CaComp", changing_pars = "B") 
    print 'starting'
    for stimulus in linspace(0.0,1.0,200):
        s += ms.solve(par_values = [stimulus])

def test2():
    s = Solutions("CICR model: Effect of stimulus on citosolic calcium")
    for stimulus in linspace(0.0,1.0,200):
        m.B = stimulus
        s += solve(m, tf = 6.0, npoints = 1000, title = 'stimulus = %g'% (m.B), outputs="Ca CaComp")#mytransformation)

def profile_test():
    # This is the main function for profiling 
    import cProfile, pstats, StringIO
    prof = cProfile.Profile()
    prof = prof.runctx("test2()", globals(), locals())
    stream = StringIO.StringIO()
    stats = pstats.Stats(prof, stream=stream)
    stats.sort_stats("time")  # Or cumulative
    stats.print_stats(40)  # 40 = how many to print
    # The rest is optional.
    #stats.print_callees()
    #stats.print_callers()
    print stream.getvalue()
    #logging.info("Profile data:\n%s", stream.getvalue())


if __name__ == "__main__":
##     profile_test()
    run_normal()

