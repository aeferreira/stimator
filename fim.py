#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
from model import *

def test():
    
    def force(A, t):
        return t*A
    m = Model('My first model')
    m.v1 = "A+B -> C", 3
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = "B   ->  "  , "2*input1"
    m.D  = variable("-kout * D")
    m.kout = 2
    m.B  = 2.2
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.Km3 = 1,6
    m.input1 = forcing(force)
    print m
        
    print '\nJacobian ----------------------------------------'
    dxdtstrings = m.dXdt_strings()
    for d in dxdtstrings:
        print d
    print
    J = m.Jacobian_strings()
    nvars = len(m.variables)
    for i in range(nvars):
        for j in range(nvars):
            print (m.variables[i].name, m.variables[j].name),
            print '\t', J[i][j]
    
    print '\nJacobian evaluation -----------------------------'
    J = m.getJacobian()
    t = 0.0
    x = array([0.0, 1.0, 2.0])
    print 'with t =', t
    print 'and x =', x
    print 'Jacobian:'
    res = J(x,t)
    print res

    print '\nMatrix dfdp -------------------------------------'
    dxdtstrings = m.dXdt_strings()
    for d in dxdtstrings:
        print d
    print
    pars = ["Km3", "V3", "kout"]
    print 'parameters are', pars
    J = m.dfdp_strings(pars)
    nvars = len(m.variables)
    for i in range(nvars):
        for j in range(len(pars)):
            print (m.variables[i].name, pars[j]),
            print '\t', J[i][j]
    

if __name__ == "__main__":
    test()

