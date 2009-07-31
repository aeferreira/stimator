#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from time import time
import de
from numpy import *
from scipy import integrate
from model import *
import timecourse

def simpleSSD(differences):
    return sum(differences*differences)
