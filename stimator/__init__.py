"""S-timator package"""

from stimator.timecourse import readTCs, read_tc, Solution, Solutions, TimeCourses

from stimator.model import Model
from stimator.modelparser import read_model


__version__ = '0.9.135'


class VersionObj(object):
    def __init__(self):
        self.version = __version__
        self.fullversion = self.version
        self.date = "Dec 2020"

    def __str__(self):
        return self.version


__full_version__ = VersionObj()


if __name__ == '__main__':
    print(__version__)
