#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator timecourse functions
# Copyright António Ferreira 2006-2009
#----------------------------------------------------------------------------
"""This module contains code to parse a model definition text,

The most important class ('StimatorParser') holds data to represent a valid model.
The parsing loop relies on regular expressions."""

import re
import math
from numpy import *

#----------------------------------------------------------------------------
#         TIME COURSE READING FUNCTION
#----------------------------------------------------------------------------

fracnumberpattern = r"[-]?\d*[.]?\d+"
realnumberpattern = fracnumberpattern + r"(e[-]?\d+)?"
identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)
realnumber = re.compile(realnumberpattern, re.IGNORECASE)

def readTimeCourseFromFile(file, atindexes=None):
    """Reads a time course from file.
    
    Returns a header with variable names (possibly empty) and a 2D numpy array with data.
    """
    
    header = []
    nvars = 0
    rows = []
    headerFound = False
    t0found = False

    f = file
    isname = False
    try:                                  
        f = open(f) # could be a name,instead of an open file
        isname = True
    except (IOError, OSError, TypeError):            
        pass                              

    for line in f:
        line = line.strip()
        if len(line) == 0:continue          #empty lines are skipped
        if line.startswith('#'): continue   #comment lines are skipped
        
        items = line.split()
        
        if identifier.match(items[0]):
            if not headerFound and not t0found:
                header = filter (identifier.match, items)
                headerFound = True
            else:
                continue
        elif not realnumber.match(items[0]):
            continue
        else:
            if not t0found:
                nvars = len(items)
                t0found = True
            temprow = [nan]*nvars
            for (i,num) in enumerate(items):
                if realnumber.match(num):
                    if atindexes:
                        temprow[atindexes[i]] = float(num)
                    else:
                        temprow[i] = float(num)
            rows.append(temprow)
    if isname:
        f.close()
    
    return header, array(rows)

class TimeCourseCollection(object):
    def __init__(self):
        self.reset()
    
    def reset(self):
        self.data = []
        self.headers = []
        self.shapes = []
        self.shortnames = []
        self.filenames = []

#----------------------------------------------------------------------------
#         TESTING CODE
#----------------------------------------------------------------------------

if __name__ == "__main__":

    print '\n===Parsing in-code timecourse ========================'

    demodata = """
#this is demo data with a header
t x y z
0       1 0         0
0.1                  0.1

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.5  - 0.5 - -
0.6 0.6 0.8 0.9

"""
    import StringIO
    aTC = StringIO.StringIO(demodata)

    h, d =  readTimeCourseFromFile(aTC)   
    print 'source:\n------------------------------------'
    print demodata
    print '------------------------------------'
    print '\nheader:'
    print h
    print '\ndata'
    print d
    aTC.seek(0) #reset StringIO
    h, d =  readTimeCourseFromFile(aTC, atindexes=(0,3,1,2))   
    print
    print '- atindexes (0,3,1,2) --------------'
    print '\nheader:'
    print h
    print '\ndata'
    print d
    
    h, d =  readTimeCourseFromFile('examples\\TSH2b.txt')   
    print '\n\n====Parsing timecourse from file =============='
    print '\n\nData from TSH2b.txt'
    print 'header:'
    print h
    print '\ndata'
    print d
    print 'dimensions are %d by %d'% d.shape

