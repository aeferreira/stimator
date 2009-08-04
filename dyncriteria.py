#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
import timecourse


def getCriteriumFunc(weights, ydata):
    if weights is None:

        def simpleSSD(d):
            return sum(d*d)
        return simpleSSD
    return None

def WSSD(differences, W):
    return sum(differences*W*differences)
