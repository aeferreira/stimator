#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model analysis classes
# Copyright Ant�nio Ferreira 2006-2009
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
    f = model.getdXdt(scale)
    #get initial values, possibly from a state in the model
    if isinstance(initial, str) or isinstance(initial, StateArray):
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
        f = model.genTransformationFunction(outputs)
        sol.apply_transf(f, f.names)
    else:
        raise TypeError("'outputs' parameter is of the wrong type")
    return sol

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
                raise ValueError( "No data for '%s' in timecourse" % str(key))
            return self.data.__getitem__(i)
        return self.data.__getitem__(key)
    def state_at(self, t):
        """Retrieves a State object with values at a time point.
        
           May have to interpolate."""
        if t > self.t[-1] or t < self.t[0]:
            raise ValueError( "No data for time '%s' in timecourse" % str(t) )
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

    def i_time(self,t):
        """Retrieves the closest index for time t."""
        if t > self.t[-1] or t < self.t[0]:
            raise ValueError( "No data for time '%s' in timecourse" % str(t) )
        # Find closest:
        ileft  = self.t.searchsorted(t, side = 'left')
        iright = self.t.searchsorted(t, side = 'right')
        if iright == ileft:
            ileft -= 1
            tl = self.t[ileft]
            tr = self.t[iright]
            if (t-tl) <= (tr-t):
                return ileft
            else:
                return iright
        else:
            return ileft
        
    def __getLastState(self):
        """Retrieves state_at last timepoint"""
        y = self.data[:,-1]
        return StateArray(dict([(x, value) for (x, value) in zip(self.names, y)]), '?')    
    last = property(__getLastState) #'last' is a synonymous, used as 'sol.last'
    def __getNumberOfTimes(self):
        """Retrieves the number of time points"""
        return self.data.shape[1]
    ntimes = property(__getNumberOfTimes)
    
    def apply_transf(self,f, newnames=None):
        """Applies a transformation to time series.
        
           f is the transformation function, with signature
           f(variables,t). variables is an array, list or tuple, t is a scalar.
           This function must return an array with the same size as 'variables'.
           newnames is a list of names of the transformed variables.
           results are kept 'in place': data is substituted."""
           
        def newf(newdata,f):
            return f(newdata[1:], newdata[0])
        trf   = apply_along_axis(newf, 0, vstack((self.t,self.data)), f)
        if newnames is not None:
            self.names = newnames
        self.data = trf

def transform(solution, f, outputs=False, title = None):
    pass
    #~ times = copy(solution.t)
    
    #~ sol = SolutionTimeCourse (times, Y.T, names, title)


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
                    raise TypeError( "Must add a solutions or collection of solutions")
            self.solutions.extend(list(other))
        elif isinstance(other, SolutionTimeCourse):
            self.solutions.append(other)
        else:
            raise TypeError( "Must add a solutions or collection of solutions")
        return self
    def __iter__(self):
        return iter(self.solutions)
    def append(self, other):
        return self.__iadd__(other)


def plot(solutions, figure = None, style = None, titles=None, ynormalize = False, superimpose = False):
    p.figure()
    colours = ['r-', 'b-', 'g-', 'k-', 'y-']
    ntc = len(solutions)
    ncols = int(math.ceil(math.sqrt(ntc)))
    nrows = int(math.ceil(float(ntc)/ncols))
    first = True
    
    if superimpose:
        p.subplot(1,1,1)
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
                p.plot(solution.t, solution[i], colours[icolour], label = label)
                icolour +=1
                if icolour == len(colours):
                    icolour = 0
        p.grid()
        p.legend(loc='best')
        p.xlabel('')
        p.ylabel('')
        p.title(solutions.title)
    else:
        for isolution,solution in enumerate(solutions):
            p.subplot(nrows,ncols,isolution+1)
            icolour = 0
            rangelines = range(len(solution))
            names = ['n/a' for i in rangelines]
            for i, name in enumerate(solution.names):
                names[i] = name
            for i in rangelines:
                p.plot(solution.t, solution[i], colours[icolour], label=names[i])
                icolour += 1
                if icolour == len(colours):
                    icolour = 0 
            p.grid()
            p.legend(loc='best')
            p.xlabel('')
            p.ylabel('')
            if titles is not None:
                p.title(titles[isolution])
            else:
                p.title(solution.title)
            yscale = p.ylim()
            if first:
                yscale_all = list(yscale)
                first = False
            else:
                if yscale[0] < yscale_all[0]: yscale_all[0] = yscale[0]
                if yscale[1] > yscale_all[1]: yscale_all[1] = yscale[1]

    if not superimpose and ynormalize:
        for isolution in range(ntc):
            p.subplot(nrows,ncols,isolution+1)
            p.ylim(yscale_all)
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

    print '--- Last time point ----'
    print 'At t =', solution1.t[-1]
    for x,value in solution1.last:
        print "%-8s= %f" % (x, value)

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

    #print '---------------- EXAMPLE 3 ------------------'
    m3 = Model("Calcium Spikes")
    m3.v0 = " -> Ca", 1
    m3.v1 = react(" -> Ca", rate = "B*k1")
    m3.k1 = 7.3
    m3.B  = 0.4
    m3.export = " Ca -> ", 10
    m3.leak   = "CaComp -> Ca", 1
    m3.v2     = "Ca -> CaComp", "65 * Ca**2 / (1+Ca**2)"
    m3.v3     = "CaComp -> Ca", "500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)"
    m3.init   = state(Ca = 0.1, CaComp = 0.63655)

    #print m3

    solution3 = solve(m3, tf = 8.0, npoints = 2000)

    print '---------------- EXAMPLE 4 ------------------'
    m4 = Model("Rossler")
    m4.X1 = variable("X2 - X3")
    m4.X2 = variable("0.36 * X2 - X1")
    m4.X3 = variable("X1 *X3 - 22.5 * X3 - 49.6 * X1 + 1117.8")
    m4.x3 = transf("X3 -50.0")
    m4.x1 = transf("X1 -18.0")
    m4.x2 = transf("X2 -50.0")
    m4.init = state(X1 = 19.0, X2 = 47, X3 = 50)

    print m4

    solution4 = solve(m4, tf = 100.0, npoints = 2000, outputs="x1 x2 x3")
    
    def transformation(vars,t):
        if t > 40.0:
            return (vars[0]-5.0, vars[1], vars[2])
        else:
            return (-5.0, vars[1], vars[2])

    solution4.apply_transf(transformation)

    plot ([solution1, solution2, solution3, solution4])

if __name__ == "__main__":
    test()