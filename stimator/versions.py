versions = {}
# populate versions dictionary

# these are standard
import sys
import platform
versions['platform'] = platform.platform()
versions['python'] = sys.version

try:
    import numpy
    versions['numpy'] = numpy.version.version
except ImportError:
    versions['numpy'] = None

try:
    import scipy
    versions['scipy'] = scipy.version.version
except ImportError:
    versions['scipy'] = None

try:
    import matplotlib
    versions['matplotlib'] = matplotlib.__version__
except ImportError:
    versions['matplotlib'] = None

try:
    import pandas
    versions['pandas'] = pandas.__version__
except ImportError:
    versions['pandas'] = None

try:
    import sympy
    versions['sympy'] = sympy.__version__
except ImportError:
    versions['sympy'] = None

try:
    import IPython
    versions['IPython'] = IPython.__version__
except ImportError:
    versions['IPython'] = None

try:
    import stimator
    versions['S-timator'] = stimator.__version__
except ImportError:
    versions['S-timator'] = None


def version_info(plist=None):
    fmtstring = '%-16s : %s\n'
    if plist is None or plist == 'all':
        plist = ['platform', 'python', 'numpy',
                 'scipy', 'matplotlib', 'pandas',
                 'IPython', 'S-timator']
    elif plist in ['packages', 'modules']:
        plist = ['numpy', 'scipy', 'matplotlib', 'pandas',
                 'IPython', 'S-timator']
    elif plist in ['platform', 'python']:
        plist = ['platform', 'python']
    else:
        pass
    
    res = ''
    for p in plist:
        v = versions[p]
        if v is None:
            res += fmtstring % (p, 'not installed')
        else:
            res += fmtstring % (p, v)
    return res

if __name__ == '__main__':
    print (version_info())
