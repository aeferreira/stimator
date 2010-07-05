#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2010 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *


def _transform2array(vect):
    if isinstance(vect, float) or isinstance(vect, int):
        res = array((vect), dtype=float)
    elif isinstance(vect, list) or isinstance(vect, tuple):
        res = diag(array(vect, dtype=float))
    else:
        res = vect # is already an array (must be 2D)
    return res
    

def constError_func(vect):
    res = _transform2array(vect)
    def CE(x):
        return res
    return CE

def propError_func(vect):
    res = _transform2array(vect)
    def CE(x):
        return float(vect) * x
    return CE
