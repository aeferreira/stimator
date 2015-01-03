"""S-timator package"""

from model import Model, register_kin_func
from dynamics import solve, ModelSolver
from timecourse import readTCs, Solutions, TimeCourses
from modelparser import read_model
import examples

class VersionObj(object):
    pass

__version__ = VersionObj()

__version__.version = '0.9.90'
__version__.fullversion = __version__.version
__version__.date = "Jun 2014"
