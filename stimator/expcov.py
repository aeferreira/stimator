#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2010 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *

def constError_func(vect):
    if isinstance(vect, float) or isinstance(vect, int): #constant for all variables
        res = array((vect), dtype=float)
    else:
        res = diag(array(vect, dtype=float))
    def CE(x):
        return res
    return CE

def propError_func(vect):
    def CE(x):
        return float(vect) * x
    return CE
