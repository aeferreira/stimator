"""S-timator package"""

from model import Model, react, transf, state, variable, register_kin_func, uncertain, varnames, variables, extvariables, transformations, reactions, parameters
from analysis import solve, plot, ModelSolver
from timecourse import readTCs, Solutions, TimeCourses
from modelparser import read_model


__version__ = '0.9.8.'+"$Revision$"[10:-1].strip()
__versiondate__ = "$Date$"[6:-1].strip()