"""S-timator package"""
from __future__ import absolute_import

from stimator.model import Model
from stimator.dynamics import solve
from stimator.timecourse import readTCs, read_tc, Solution, Solutions, TimeCourses
from stimator.modelparser import read_model
import stimator.examples as examples

class VersionObj(object):
    pass

__version__ = VersionObj()

__version__.version = '0.9.99'
__version__.fullversion = __version__.version
__version__.date = "Mar 2016"
