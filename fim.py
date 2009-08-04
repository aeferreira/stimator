#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
from model import *
import re
import sympy

identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)

def identifiersInExpr(_expr):
    iterator = identifier.finditer(_expr)
    return [_expr[m.span()[0]:m.span()[1]] for m in iterator]

def dexpr_dname(_expr, _name):
    _ids = identifiersInExpr(_expr)
    if _name not in _ids:
        return '0'
    


def test():
    
    def force(A, t):
        return t*A
    m = Model('My first model')
    m.v1 = "A+B -> C", 3
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = "B   ->  "  , "2*input1"
    m.D  = variable("-2 * D")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.Km3 = 1,6
    m.init = state(A = 1.0, C = 1, D = 1)
    m.afterwards = state(A = 1.0, C = 2, D = 1)
    m.afterwards.C.uncertainty(1,3)
    m.input1 = forcing(force)
    

    print '********** Testing rate and dXdt strings *******************'
    print 'rates_strings(): -------------------------'
    for v in m.rates_strings():
        print v
    print '\nids in rates -----------------------------'
    for v in m.rates_strings():
        ids = identifiersInExpr(v[1])
        print v[0], ids

if __name__ == "__main__":
    test()

