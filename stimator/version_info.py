"""S-timator package"""
class VersionObj(object):
    pass

__version__ = VersionObj()

__version__.version = '0.9.85'
__version__.fullversion = __version__.version+"$Revision$"[10:-1].strip()
__version__.date = "$Date$"[6:-1].strip()

