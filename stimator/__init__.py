"""S-timator package"""

from model import Model
from dynamics import solve
from timecourse import readTCs, read_tc, Solution, Solutions, TimeCourses
from modelparser import read_model
import examples

class VersionObj(object):
    pass

__version__ = VersionObj()

__version__.version = '0.9.99'
__version__.fullversion = __version__.version
__version__.date = "Mar 2016"
