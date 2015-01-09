"""S-timator : DEMO of scanning parameter functions."""
from stimator import *
from stimator.dynamics import ModelSolver
from time import time, sleep
from numpy import append, linspace

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

def run_normal():
    npoints = 1000
    npars = 200
    title = "CICR model: Effect of stimulus on citosolic calcium"
    print (title)
    print 'stressing solve() and ModelSolver.solve()'
    print 'with %d time points and %d parameters' % (npoints, npars)
    s = Solutions(title=title)
    
    time0 = time()
    ms = ModelSolver(m,tf=6.0, npoints=npoints, 
                     outputs="Ca CaComp", 
                     changing_pars="B") 
    print '\nstarting'
    for stimulus in linspace(0.0,1.0,npars):
        s += ms.solve(title='stimulus = %g' % stimulus, par_values=[stimulus])

    print '\nusing ModelSolver done in', time()-time0, 's'

    s = Solutions(title=title)

    time0 = time()
    for stimulus in linspace(0.0,1.0,npars):
        m.parameters.B = stimulus
        s += m.solve(tf=6.0, npoints=npoints, 
                    title='stimulus = %g'% (m.parameters.B), 
                    outputs="Ca CaComp")#mytransformation)

    print 'using solve done in', time()-time0, 's'

    s = Solutions(title="CICR model: compare ModelSolver and solve()")
    svalues = (0.2, 0.4, 0.78)

    ms = ModelSolver(m,tf=6.0, npoints=npoints, 
                     outputs="Ca CaComp", changing_pars = "B") 
    for stimulus in svalues:
        s += ms.solve(title='stimulus = %g' % stimulus, par_values=[stimulus])

    for stimulus in svalues:
        m.parameters.B = stimulus
        s += m.solve(tf=6.0, npoints=npoints, 
                    title='stimulus = %g'% (m.parameters.B), 
                    outputs="Ca CaComp")#mytransformation)

    s.plot(ynormalize = True, fig_size=(16,9), show = True)
    #plot(s, superimpose=True)

def test():
    s = Solutions("CICR model: Effect of stimulus on citosolic calcium")
    ms = ModelSolver(m,tf=6.0, npoints=1000, 
                     outputs="Ca CaComp", changing_pars = "B") 
    print 'starting'
    for stimulus in linspace(0.0,1.0,200):
        s += ms.solve(par_values=[stimulus])

def test2():
    s = Solutions("CICR model: Effect of stimulus on citosolic calcium")
    for stimulus in linspace(0.0,1.0,200):
        m.parameters.B = stimulus
        s += m.solve(tf=6.0, npoints=1000, 
                    title='stimulus = %g'% (m.parameters.B), 
                    outputs="Ca CaComp")#mytransformation)

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

