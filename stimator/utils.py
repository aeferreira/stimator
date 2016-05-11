"""util functions:

   - tests for common types
   - function argument parsing
   - simple format of h:m:s from seconds
"""

from __future__ import print_function, absolute_import
from six import string_types, integer_types


def _args_2_dict(*p, **pdict):
    """Transform arguments to a dict, as in dict() plus f(a,b) -> {a:b}."""
    if len(p) == 2:
        p = ({p[0]: p[1]},)
    dpars = dict(*p, **pdict)
    return dpars


def make_dict_from_args(*args, **kwargs):
    return dict(*args, **kwargs)

def _is_sequence(arg):
    isstring = isinstance(arg, string_types)
    isothersequences = hasattr(arg, "__getitem__") or hasattr(arg, "__iter__")
    return not isstring and isothersequences

def _is_string(a):
    return isinstance(a, string_types)

def _is_number(a):
    return isinstance(a, float) or isinstance(a, integer_types)

# helper to transform string arguments in lists:
def listify(arguments):
    if isinstance(arguments, list) or isinstance(arguments, tuple):
        return [a.strip() for a in arguments]
    if _is_string(arguments):
        arguments = arguments.split()
        return [a.strip() for a in arguments]

def write2file(filename, astring):
    f = open(filename, 'w')
    f.write(astring)
    f.close()

def readfile(filename):
    f = open(filename)
    s = f.read()
    f.close()
    return s

def s2HMS(seconds):
    m, s = divmod(seconds, 60.0)
    h, m = divmod(m, 60.0)
    if h == 0:
        return "%02dm %06.3fs" % (m, s)
    return "%dh %02dm %06.3fs" % (h, m, s)

if __name__ == '__main__':
    for t in ['ok ok', u'pppp', 123.4, 12, (1, 2, 'oo'), [4, 5, 'OO']]:
        if _is_string(t):
            print('{} is a string'.format(t))
        else:
            print('{} is not a string'.format(t))

        if _is_number(t):
            print('{} is a number'.format(t))
        else:
            print('{} is not a number'.format(t))

        if _is_sequence(t):
            print('{} is a sequence'.format(t))
        else:
            print('{} is not a sequence'.format(t))
        print()
