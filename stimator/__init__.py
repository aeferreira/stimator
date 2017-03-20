"""S-timator package"""
from __future__ import print_function, absolute_import

from stimator.model import Model
from stimator.dynamics import solve
from stimator.timecourse import readTCs, read_tc, Solution, Solutions, TimeCourses
from stimator.modelparser import read_model
import stimator.examples as examples

class VersionObj(object):
    def __init__(self):
        self.version = '0.9.110'
        self.fullversion = self.version
        self.date = "Mar 2017"

    def __str__(self):
        return self.version

__version__ = VersionObj()


if __name__ == '__main__':
    print(__version__)