"""S-timator package"""
from __future__ import print_function, absolute_import
import os
import sys

from stimator.timecourse import readTCs, read_tc, Solution, Solutions, TimeCourses
#from stimator.dynamics import solve
from stimator.model import Model
from stimator.modelparser import read_model
import stimator.examples as examples

__version__ = '0.9.120'

sys.path.insert(0, os.path.split(os.path.split(os.path.abspath(__file__))[0])[0])

from tests import run_tests

class VersionObj(object):
    def __init__(self):
        self.version = __version__
        self.fullversion = self.version
        self.date = "Oct 2018"

    def __str__(self):
        return self.version

__full_version__ = VersionObj()


if __name__ == '__main__':
    print(__version__)