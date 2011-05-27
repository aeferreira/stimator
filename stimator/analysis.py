#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model analysis classes
# Copyright António Ferreira 2006-2009
#----------------------------------------------------------------------------
import math
from model import *
from dynamics import *
from modelparser import read_model
from numpy import *
from scipy import integrate
from timecourse import SolutionTimeCourse, Solutions

import pylab as p

def solve(model, tf = 1.0, npoints = 500, t0 = 0.0, initial = 'init', times = None, outputs=False, title = None):
    salg=integrate._odepack.odeint
    names = [x for x in varnames(model)]

    #get initial values, possibly from a state in the model
    if isinstance(initial, str) or isinstance(initial, StateArray):
        y0 = copy(state2array(model,initial))
    else:
        y0 = copy(initial)
    if times is None:
        times = linspace (t0, tf, npoints)
    # scale times to maximum time in data
    t0 = times[0]
    scale = float(times[-1] - times[0])
    #scale = 1.0
    
    f = getdXdt(model, scale=scale, t0=t0)
    t  = (times-t0)/scale  # this scales time points
    output = salg(f, y0, t, (), None, 0, -1, -1, 0, None, 
                    None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
    if output[-1] < 0: return None
    Y = output[0]
    if title is None:
        title = model.getData('title')
    
    sol = SolutionTimeCourse (times, Y.T, names, title)
    
    #get outputs
    if outputs is False: # variables are output
        pass
    elif outputs is True: #transformations are output
        #compute model transformations
        f     = model.transf_func()
        names = [x.name for x in model.transf]
        sol.apply_transf(f,names)
    elif isinstance(outputs, str) or callable(outputs): 
        #a filter string or transformation function
        f = genTransformationFunction(model, outputs)
        sol.apply_transf(f, f.names)
    else:
        raise TypeError("'outputs' parameter is of the wrong type")
    return sol

def plot(solutions, show = False, figure = None, style = None, titles=None, ynormalize = False, superimpose = False, legend=True):
    if isinstance(solutions, SolutionTimeCourse):
        s = Solutions()
        s.append(solutions)
        solutions = s
    if figure is None:
        figure = p.figure()
    colours = ['r-', 'b-', 'g-', 'k-', 'y-']
    ntc = len(solutions)
    ncols = int(math.ceil(math.sqrt(ntc)))
    nrows = int(math.ceil(float(ntc)/ncols))
    first = True
    
    if superimpose:
        curraxis=figure.add_subplot(1,1,1)
        icolour = 0                
        for isolution,solution in enumerate(solutions):
            rangelines = range(len(solution))
            names = ['n/a' for i in rangelines]
            for i, name in enumerate(solution.names):
                names[i] = name
            for i in rangelines:
                if len(solution) == 1:
                    label = "%s"%(solution.title)
                else:
                    label = "%s, %s"%(names[i], solution.title)
                curraxis.plot(solution.t, solution[i], colours[icolour], label = label)
                icolour +=1
                if icolour == len(colours):
                    icolour = 0
        curraxis.grid()
        h, l = curraxis.get_legend_handles_labels()
        curraxis.legend(h, l, loc='best')
        curraxis.set_xlabel('')
        curraxis.set_ylabel('')
        if hasattr(solutions, 'title'):
            curraxis.set_title(solutions.title)
    else:
        for isolution,solution in enumerate(solutions):
            curraxis=figure.add_subplot(nrows,ncols,isolution+1)
            icolour = 0
            rangelines = range(len(solution))
            names = ['n/a' for i in rangelines]
            for i, name in enumerate(solution.names):
                names[i] = name
            for i in rangelines:
                curraxis.plot(solution.t, solution[i], colours[icolour], label=names[i])
                icolour += 1
                if icolour == len(colours):
                    icolour = 0 
            curraxis.grid()
            h, l = curraxis.get_legend_handles_labels()
            curraxis.legend(h, l, loc='best')

            curraxis.set_xlabel('')
            curraxis.set_ylabel('')
            if titles is not None:
                curraxis.set_title(titles[isolution])
            else:
                curraxis.set_title(solution.title)
            yscale = curraxis.get_ylim()
            if first:
                yscale_all = list(yscale)
                first = False
            else:
                if yscale[0] < yscale_all[0]: yscale_all[0] = yscale[0]
                if yscale[1] > yscale_all[1]: yscale_all[1] = yscale[1]

    if not superimpose and ynormalize:
        for isolution in range(ntc):
            curraxis=figure.add_subplot(nrows,ncols,isolution+1)
            curraxis.set_ylim(yscale_all)
    if show:
        p.show()

def test():
    print '---------------- EXAMPLE 1 ------------------'
    m1 = read_model("""
    title Glyoxalase system
    glo1 = HTA -> SDLTSH, rate = V1*HTA/(Km1 + HTA)
    glo2 = SDLTSH ->    , rate = V2*SDLTSH/(Km2 + SDLTSH)
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
    v1 = A -> B, rate = k1*A
    k1 = 10
    v2 = B -> C, rate = k2*B**0.5
    k2 = 5
    v3 = C -> D, rate = k3*C**0.5
    k3 = 2
    v4 = C -> E, rate = k4*C**0.5
    k4 = 8
    v5 = D ->  , rate = k5*D**0.5
    k5 = 1.25
    v6 = E ->  , rate = k6*E**0.5
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

    plot ([solution1, solution2, solution3, solution4], superimpose=False, show = True)

if __name__ == "__main__":
    test()
