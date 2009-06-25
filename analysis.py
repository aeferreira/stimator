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
from numpy import *
from scipy import integrate

import pylab as p

def solve(model, tf = 1.0, npoints = 500, t0 = 0.0, initial = 'init', times = None, outputs=False, title = None):
    salg=integrate._odepack.odeint
    names = [x.name for x in model.variables]
    scale = 1.0
    f = model.scaled_dXdt(scale)
    if isinstance(initial, str) or isinstance(initial, State):
        y0 = copy(model.vectorize(initial))
    else:
        y0 = copy(initial)
    if times is None:
        times = linspace (t0, tf, npoints)
    t  = copy(times)
    output = salg(f, y0, t, (), None, 0, -1, -1, 0, None, 
                    None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
    if output[-1] < 0: return None
    Y = output[0]
    if title is None:
        title = model.title
    sol = SolutionTimeCourse (times, Y.T, names, title)
    if outputs is False:
        return sol
    elif outputs is True:
        #compute transf
        f     = model.transf_func()
        def newf(newdata,f):
            return f(newdata[1:], newdata[0])
        names = [x.name for x in model.transf]
        trf   = apply_along_axis(newf, 0, vstack((sol.t,sol.data)), f)
        return SolutionTimeCourse (sol.t, trf, names, title)
    elif isinstance(outputs, str) or callable(outputs):
        f = model.genTransformationFunction(outputs)
        def newf(newdata,f):
            return f(newdata[1:], newdata[0])
        trf   = apply_along_axis(newf, 0, vstack((sol.t,sol.data)), f)
        if not callable(outputs):
            names = outputs.split()
        else:
            if hasattr(outputs, 'names'):
                names = outputs.names
            else:
                names = ['transf%d'% i for i in range(trf.shape[0])]
        return SolutionTimeCourse (sol.t, trf, names, title)

class SolutionTimeCourse(object):
    """Holds a timecourse created by ODE solvers"""
    def __init__(self, t = array([]), data = array([]), names = [], title = ""):
        self.t = t         #values of time points
        self.data = data   # table of solution points: series in rows, times in cols
        self.names = names # names of the series
        self.shape = data.shape
        self.title = title # a title for the solution
        
    def __len__(self):
        return self.data.shape[0]
    def __nonzero__(self):
        return len(t) > 0
    def __getitem__(self, key):
        """retrieves a series by name or index"""
        if isinstance(key, str) or isinstance(key, unicode):
            try:
                i = self.names.index(key)
            except ValueError:
                raise ValueError, "No data for '%s' in timecourse" % str(key)
            return self.data.__getitem__(i)
        return self.data.__getitem__(key)
    def state_at(self, t):
        """Retrieves a State object with values at a time point.
        
           May have to interpolate."""
        if t > self.t[-1] or t < self.t[0]:
            raise ValueError, "No data for time '%s' in timecourse" % str(t)
        # Interpolate:
        ileft = self.t.searchsorted(t, side = 'left')
        iright = self.t.searchsorted(t, side = 'right')
        if iright == ileft:
            ileft -= 1
            tl = self.t[ileft]
            tr = self.t[iright]
            yl = self.data[:,ileft]
            yr = self.data[:,iright]
            m = (yr-yl)/(tr-tl)
            y = yl + m *(t-tl)
        else:
            y = self.data[:, ileft]
        return StateArray(dict([(x, value) for (x, value) in zip(self.names, y)]), '?')
    def __getLastState(self):
        """Retrieves state_at last timepoint"""
        y = self.data[:,-1]
        return StateArray(dict([(x, value) for (x, value) in zip(self.names, y)]), '?')    
    last = property(__getLastState) #'last' is a synonymous, used as 'sol.last'

    
class Solutions(object):
    """Holds a colection of objects of class SolutionTimeCourse"""
    def __init__(self, title = ""):
        self.title = title
        self.solutions = []
    
    def __getitem__(self, key):
        """retrieves a series by index"""
        return self.solutions.__getitem__(key)
    def __len__(self):
        return len(self.solutions)
    def __nonzero__(self):
        return len(self.solutions) > 0
    def __iadd__(self,other):
        if isinstance(other, Solutions):
            self.solutions.extend(other.solutions)
        elif isinstance(other, list) or isinstance(other, tuple):
            for s in other:
                if not isinstance(s, SolutionTimeCourse): 
                    raise TypeError, "Must add a solutions or collection of solutions"
            self.solutions.extend(list(other))
        elif isinstance(other, SolutionTimeCourse):
            self.solutions.append(other)
        else:
            raise TypeError, "Must add a solutions or collection of solutions"
        return self
    def __iter__(self):
        return iter(self.solutions)
    def append(self, other):
        return self.__iadd__(other)


def plot(solutions, figure = None, style = None, titles=None):
    p.figure()
    colours = ['r-', 'b-', 'g-', 'k-', 'y-']
    ntc = len(solutions)
    ncols = int(math.ceil(math.sqrt(ntc)))
    nrows = int(math.ceil(float(ntc)/ncols))

    for isolution,solution in enumerate(solutions):
        p.subplot(nrows,ncols,isolution+1)
        for i in range(len(solution)):
            p.plot(solution.t, solution[i], colours[i], label=solution.names[i])
        p.grid()
        p.legend(loc='best')
        p.xlabel('')
        p.ylabel('')
        if titles is not None:
            p.title(titles[isolution])
        else:
            p.title(solution.title)
    p.show()

def test():
    print '---------------- EXAMPLE 1 ------------------'
    m1 = Model("Glyoxalase system")
    m1.glo1 = react("HTA -> SDLTSH", rate = "V1*HTA/(Km1 + HTA)")
    m1.glo2 = react("SDLTSH -> "   , rate = "V2*SDLTSH/(Km2 + SDLTSH)")
    m1.V1  = 2.57594e-05
    m1.Km1 = 0.252531
    m1.V2  = 2.23416e-05
    m1.Km2 = 0.0980973
    m1.init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)

    print m1

    solution1 = solve(m1, tf = 4030.0)

    #print t
    #print solution
    #print trf

    print '--- Last time point ----'
    print 'At t =', solution1.t[-1]
    for x,value in solution1.last:
        print "%-8s= %f" % (x, value)

    #plot results...
    #~ f1 = p.figure(1)
    #~ p.subplot(221) 
    #~ plot(solution, title = m1.title)

    #print '---------------- EXAMPLE 2 ------------------'
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

    #print m2
    times = append(linspace(0.0,5.0,500),linspace(5.0,10.0, 500))

    solution2 = solve(m2, tf = 10.0, times=times)
    #~ p.subplot(222) 
    #~ plot(solution, title = m2.title)

    #print '---------------- EXAMPLE 3 ------------------'
    m3 = Model("Calcium Spikes")
    m3.v0 = react(" -> Ca", 1)
    m3.v1 = react(" -> Ca", rate = "B*k1")
    m3.k1 = 7.3
    m3.B  = 0.4
    m3.export = react(" Ca -> ", 10)
    m3.leak   = react("CaComp -> Ca", 1)
    m3.v2     = react("Ca -> CaComp", rate = "65 * Ca**2 / (1+Ca**2)")
    m3.v3     = react("CaComp -> Ca", rate = "500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)")
    m3.init   = state(Ca = 0.1, CaComp = 0.63655)

    #print m3

    solution3 = solve(m3, tf = 8.0, npoints = 2000)
    #~ p.subplot(223) 
    #~ plot(solution, title = m3.title)

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

    solution4 = solve(m4, tf = 100.0, npoints = 2000, outputs="x1 x2 x3")

    plot ([solution1, solution2, solution3, solution4])

if __name__ == "__main__":
    test()
