#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator timecourse functions
# Copyright António Ferreira 2006-2009
#----------------------------------------------------------------------------
import os.path
import re
import math
from numpy import *
import model

#----------------------------------------------------------------------------
#         TIME COURSE READING FUNCTION
#----------------------------------------------------------------------------

fracnumberpattern = r"[-]?\d*[.]?\d+"
realnumberpattern = fracnumberpattern + r"(e[-]?\d+)?"
identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)
realnumber = re.compile(realnumberpattern, re.IGNORECASE)

def readTimeCourseFromFile(file, atindexes=None):
    """Reads a time course from file.
    
    Returns a header with variable names (possibly absent in file) and a 2D numpy array with data.
    """
    
    header = []
    nvars = 0
    rows = []
    headerFound = False
    t0found = False

    f = file
    isname = False
    try:                                  
        f = open(f, "rU") # could be a name,instead of an open file
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
    rows = array(rows)
    #create default names "t, x1, x2, x3,..."
    if len(header) == 0:
        header = ['t']
        for i in range(1, rows.shape[1]):
            header.append('x%d'%i)

    return header, rows



class TimeCourseCollection(object):
    def __init__(self):
        self.reset()
    
    def reset(self):
        self.data           = []
        self.headers        = []
        self.shapes         = []
        self.shortnames     = []
        self.filenames      = []
        self.basedir        = None
        self.intvarsorder   = None
        self.variablesorder = None # list of names indicating the order of variables in timecourses

    def loadTimeCourses (self,filedir):

        if len(self.filenames) == 0 :
           print "No time courses to load!\nPlease indicate some time courses with 'timecourse <filename>'"
           return 0
        
        # check and load timecourses
        self.basedir = filedir
        cwd = os.getcwdu()
        os.chdir(self.basedir)
        pathlist = [os.path.abspath(k) for k in self.filenames]

        print "-------------------------------------------------------"
        self.data = []
        for filename in pathlist:
            if not os.path.exists(filename) or not os.path.isfile(filename):
                print "Time course file \n%s\ndoes not exist"% filename
                os.chdir(cwd)
                return 0
            h,d = readTimeCourseFromFile(filename, atindexes=self.intvarsorder)
            if d.shape == (0,0):
                print "File\n%s\ndoes not contain valid time-course data"% filename
                os.chdir(cwd)
                return 0
            else:
                print "%d time points for %d variables read from file %s" % (d.shape[0], d.shape[1]-1, filename)
                self.headers.append(h)
                self.data.append(d)
        #~ for i,d in enumerate(self.data):
            #~ if d.shape[1] != len(self.model.variables)+1:
                #~ print "There are %i initial values in time course %s but model has %i variables"%(d.shape[1]-1,
                               #~ self.tc.filenames[i],len(self.model.variables))
                #~ return None
        self.shapes     = [i.shape for i in self.data]
        self.shortnames = [os.path.split(filename)[1] for filename in pathlist]
        os.chdir(cwd)
        return len(pathlist)

def readTimeCourses(filenames, filedir, intvarsorder):
    tc = TimeCourseCollection()
    tc.filenames = filenames
    tc.intvarsorder = intvarsorder
    nread = tc.loadTimeCourses(filedir)
    return tc


class SolutionTimeCourse(object):
    """Holds a timecourse created by ODE solvers"""
    def __init__(self, t = array([]), data = array([]), names = [], title = ""):
        self.t = t         #values of time points
        self.data = data   # table of solution points: series in rows, times in cols
        self.names = names # names of the series
        self.shape = data.shape
        self.title = title # a title for the solution
        
    def __len__(self):
        return self.data.shape[0]
    def __nonzero__(self):
        return len(t) > 0
    def __getitem__(self, key):
        """retrieves a series by name or index"""
        if isinstance(key, str) or isinstance(key, unicode):
            try:
                i = self.names.index(key)
            except ValueError:
                raise ValueError( "No data for '%s' in timecourse" % str(key))
            return self.data.__getitem__(i)
        return self.data.__getitem__(key)
    def state_at(self, t):
        """Retrieves a State object with values at a time point.
        
           May have to interpolate."""
        if t > self.t[-1] or t < self.t[0]:
            raise ValueError( "No data for time '%s' in timecourse" % str(t) )
        # Interpolate:
        ileft = self.t.searchsorted(t, side = 'left')
        iright = self.t.searchsorted(t, side = 'right')
        if iright == ileft:
            ileft -= 1
            tl = self.t[ileft]
            tr = self.t[iright]
            yl = self.data[:,ileft]
            yr = self.data[:,iright]
            m = (yr-yl)/(tr-tl)
            y = yl + m *(t-tl)
        else:
            y = self.data[:, ileft]
        return model.StateArray(dict([(x, value) for (x, value) in zip(self.names, y)]), '?')

    def i_time(self,t):
        """Retrieves the closest index for time t."""
        if t > self.t[-1] or t < self.t[0]:
            raise ValueError( "No data for time '%s' in timecourse" % str(t) )
        # Find closest:
        ileft  = self.t.searchsorted(t, side = 'left')
        iright = self.t.searchsorted(t, side = 'right')
        if iright == ileft:
            ileft -= 1
            tl = self.t[ileft]
            tr = self.t[iright]
            if (t-tl) <= (tr-t):
                return ileft
            else:
                return iright
        else:
            return ileft
        
    def __getLastState(self):
        """Retrieves state_at last timepoint"""
        y = self.data[:,-1]
        return model.StateArray(dict([(x, value) for (x, value) in zip(self.names, y)]), '?')    
    last = property(__getLastState) #'last' is a synonymous, used as 'sol.last'
    def __getNumberOfTimes(self):
        """Retrieves the number of time points"""
        return self.data.shape[1]
    ntimes = property(__getNumberOfTimes)
    
    def apply_transf(self,f, newnames=None):
        """Applies a transformation to time series.
        
           f is the transformation function, with signature
           f(variables,t). variables is an array, list or tuple, t is a scalar.
           This function must return an array with the same size as 'variables'.
           newnames is a list of names of the transformed variables.
           results are kept 'in place': data is substituted."""
           
        def newf(newdata,f):
            return f(newdata[1:], newdata[0])
        trf   = apply_along_axis(newf, 0, vstack((self.t,self.data)), f)
        if newnames is not None:
            self.names = newnames
        self.data = trf
        return self

class Solutions(object):
    """Holds a colection of objects of class SolutionTimeCourse"""
    def __init__(self, title = ""):
        self.title = title
        self.solutions = []
    
    def __getitem__(self, key):
        """retrieves a series by index"""
        return self.solutions.__getitem__(key)
    def __len__(self):
        return len(self.solutions)
    def __nonzero__(self):
        return len(self.solutions) > 0
    def __iadd__(self,other):
        if isinstance(other, Solutions):
            self.solutions.extend(other.solutions)
        elif isinstance(other, list) or isinstance(other, tuple):
            for s in other:
                if not isinstance(s, SolutionTimeCourse): 
                    raise TypeError( "Must add a solutions or collection of solutions")
            self.solutions.extend(list(other))
        elif isinstance(other, SolutionTimeCourse):
            self.solutions.append(other)
        else:
            raise TypeError( "Must add a solutions or collection of solutions")
        return self
    def __iter__(self):
        return iter(self.solutions)
    def append(self, other):
        return self.__iadd__(other)


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
    print
    
    print '--Using a SolutionTimeCourse interface---------'
    sol = SolutionTimeCourse (d[:,0].T, d[:,1:].T, h[1:])
    print 'retrieving components...'
    try:
        print 'len(sol)'
        print len(sol)
        print 'sol.ntimes'
        print sol.ntimes
        print 'sol[0] (first var, "x")'
        print sol[0]
        print 'sol.t'
        print sol.t
        print "sol['x']"
        print sol['x']
        print "sol.names"
        print sol.names
        print 'Last time point, sol[:,-1] returns array'
        print sol[:,-1]
        print 'The following return model.StateArray objects:'
        print 'sol.state_at(0.2)'
        print sol.state_at(0.2)
        print 'sol.state_at(0.55)'
        print sol.state_at(0.55)
        print 'sol.state_at(0.0)'
        print sol.state_at(0.0)
        print 'sol.state_at(0.6)'
        print sol.state_at(0.6)
        print 'sol.last (Last time point the easy way)'
        print sol.last
        print 'sol.last.x'
        print sol.last.x
        print "sol['k']"
        print sol['k']
    except ValueError, msg:
        print msg
    print
    
    aTC.seek(0) #reset StringIO
    h, d =  readTimeCourseFromFile(aTC, atindexes=(0,3,1,2))   
    print
    print '- atindexes (0,3,1,2) --------------'
    print '\nheader:'
    print h
    print '\ndata'
    print d
    
    h, d =  readTimeCourseFromFile('../models/TSH2b.txt')   
    print '\n\n====Parsing timecourse from file =============='
    print '\n\nData from TSH2b.txt'
    print 'header:'
    print h
    print '\ndata'
    print d
    print 'dimensions are %d by %d'% d.shape
    
    #~ h, d =  readTimeCourseFromFile(''../models/TSH2a.txt')   
    #~ print '\n\n====Parsing timecourse from file =============='
    #~ print '\n\nData from TSH2a.txt'
    #~ print 'header:'
    #~ print h
    #~ print '\ndata'
    #~ print d
    #~ print 'dimensions are %d by %d'% d.shape

