#!/usr/bin/env python
# -*- coding: latin1-*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator timecourse functions
# Copyright António Ferreira 2006-2013
#----------------------------------------------------------------------------
import os.path
import StringIO
import re
from numpy import *
import model
import modelparser

fracnumberpattern = r"[-]?\d*[.]?\d+"
realnumberpattern = fracnumberpattern + r"(e[-]?\d+)?"
identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)
realnumber = re.compile(realnumberpattern, re.IGNORECASE)

class StimatorTCError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg


#----------------------------------------------------------------------------
#         THE BASIC TIMECOURSE CLASS
#----------------------------------------------------------------------------

class SolutionTimeCourse(object):
    """Holds a timecourse created by ODE solvers"""
    def __init__(self, t = array([]), data = array([]), names = [], title = ""):
        self.t = t         #values of time points
        self.data = data   # table of solution points: series in rows, times in cols
        self.names = names # names of the series
        self.title = title # a title for the solution
        
        #for solutions read from a file
        self.filename =  ""
        self.shortname = ""
        
        
    def __len__(self):
        """retrieves the number of vars in this solution, NOT the len(timepoints)."""
        return self.data.shape[0]
    def __nonzero__(self):
        return len(t) > 0
    
    def __getNumberOfTimes(self):
        """Retrieves the number of time points"""
        return self.data.shape[1]
    ntimes = property(__getNumberOfTimes)
    
    def __getShape(self):
        return self.data.shape
    shape = property(__getShape)
    
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
    
    def load_from_str(self, s, names = None):
        aTC   = StringIO.StringIO(s)
        aTC.seek(0)
        self.load_from(aTC, names)
    
    def load_from(self, filename, names = None):
        """Reads a time course from file.
        
        Fills self.names from a header with variable names (possibly absent in file) 
        Fills a 2D numpy array with whitespace separated data. 
        """
        
        header = []
        nvars = 0
        rows = []
        headerFound = False
        t0found = False

        if hasattr(filename, 'read'):
            f = filename
            isname = False
        else:
            f = open(filename, "rU") # could be a name,instead of an open file
            isname = True
            
        for line in f:
            line = line.strip()
            if len(line) == 0:continue          #empty lines are skipped
            if line.startswith('#'): continue   #comment lines are skipped
            #print line
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
                        temprow[i] = float(num)
                rows.append(temprow)
        if isname:
            f.close()
        
        #create default names "t, x1, x2, x3,..." or use names if provided
        if len(header) == 0:
            header = ['t']
            for i in range(1, nvars):
                header.append('x%d'%i)
            if names is not None:
                smallindx = min(len(header)-1, len(names))
                for i in range(smallindx):
                    header[i+1] = names[i]
        data = array(rows)
        self.names = header[1:]
        self.t = data[:,0].T
        self.data = data[:,1:].T

    def save_to_str(self):
        aTC   = StringIO.StringIO()
        aTC.seek(0)
        self.write_to(aTC)
        return aTC.getvalue()

    def write_to(self, filename):
        """Writes a time course to a file or file-like object.
        """
        
        header = []
        nvars = 0
        rows = []
        headerFound = False
        t0found = False

        if hasattr(filename, 'read'):
            f = filename
            isname = False
        else:
            f = open(filename, "w") # could be a name,instead of an open file
            isname = True
            
        f.write("%s %s\n"%('t', " ".join(self.names)))
        npoints = len(self.t)
        for i in range(npoints):
            row = [self.t[i]]
            row.extend(self.data[:,i])
            row = " ".join([str(j) for j in row])
            f.write("%s\n"%row)
        if isname:
            f.close()

    def clone(self):
        """Clones the entire solution."""
        tc = SolutionTimeCourse(self.t.copy(), self.data.copy(), self.names[:], self.title)
        tc.filename = self.filename
        tc.shortname = self.shortname
        return tc

    def copy(self, names = [], newtitle = None):
        """Constructs new solution, restricted to the variables in 'names'."""
        if not (isinstance(names, list) or isinstance(names, tuple)):
            names = names.strip()
            names = names.split()
        t = self.t.copy()
        if names == []:
            names = self.names
        nameindexes = []
        for name in names:
            if not name in self.names:
                raise ValueError( "No data for '%s' in timecourse" % name)
            nameindexes.append(self.names.index(name))
        data = self.data[nameindexes, :].copy()
        if newtitle is not None:
            title = newtitle
        else:
            title = self.title
        tc = SolutionTimeCourse(t, data, names[:], title)
        tc.filename = self.filename
        tc.shortname = self.shortname
        return tc
            
    def orderByNames(self, varnames):
        oldindexes = range(len(self))
        newindexes = []
##         for varname in varnames:
##             if varname not in self.names:
##                 raise StimatorTCError("series %s was not found in timecourse %s"%(varname, self.title))
        for vname in varnames:
            if vname in self.names:
                indx = self.names.index(vname)
                newindexes.append(indx)
                oldindexes.remove(indx)
        newindexes.extend(oldindexes)
        self.names = [self.names[i] for i in newindexes]
        self.data = self.data[array(newindexes, dtype=int)]

#----------------------------------------------------------------------------
#         A CONTAINER FOR TIMECOURSES
#----------------------------------------------------------------------------
class Solutions(object):
    """Holds a colection of objects of class SolutionTimeCourse"""
    def __init__(self, title = ""):
        self.title = title
        self.solutions = []
        self.shortnames     = []
        self.filenames      = []
        self.basedir        = None
        self.defaultnames   = None # list of names to use if headers are missing

    def __str__(self):
        if len(self.filenames) >0:
            return str(self.filenames)
        else:
            return 'No timecourses'
    
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
    def loadTimeCourses (self,filedir = None, names = None, verbose = False):
        if len(self.filenames) == 0 :
           print "No time courses to load!\nPlease indicate some time courses with 'timecourse <filename>'"
           return 0
        
        # check and load timecourses
        cwd = os.getcwdu()
        if filedir is not None:
            self.basedir = filedir
        else:
            self.basedir = cwd
        os.chdir(self.basedir)
        pathlist = [os.path.abspath(k) for k in self.filenames]

        self.data = []
        nTCsOK = 0
        if verbose:
            print "-------------------------------------------------------"
        for filename in pathlist:
            if not os.path.exists(filename) or not os.path.isfile(filename):
                print "Time course file \n%s\ndoes not exist"% filename
                os.chdir(cwd)
                return nTCsOK
            sol = SolutionTimeCourse()
            sol.load_from(filename, names=names)
            if sol.shape == (0,0):
                print "File\n%s\ndoes not contain valid time-course data"% filename
                os.chdir(cwd)
                return nTCsOK
            else:
                if verbose:
                    print "%d time points for %d variables read from file %s" % (sol.ntimes, len(sol), filename)
                self.append(sol)
                nTCsOK += 1
        self.shortnames = [os.path.split(filename)[1] for filename in pathlist]
        for i,sol in enumerate(self.solutions):
            sol.title = self.shortnames[i]
            sol.shortname = self.shortnames[i]
            sol.filename = self.filenames[i]
        os.chdir(cwd)
        return nTCsOK

    def saveTimeCoursesTo (self, filenames, filedir = None, verbose = False):
        if len(self) == 0 :
           print "No time courses to save!"
           return 0
                
        # check and load timecourses
        cwd = os.getcwdu()
        if filedir is not None:
            self.basedir = filedir
        else:
            self.basedir = cwd
        os.chdir(self.basedir)
        pathlist = [os.path.abspath(k) for k in filenames]

        if verbose:
            print "-------------------------------------------------------"
        for fn, sol in zip(pathlist, self.solutions):
            sol.write_to(fn)
            if verbose:
                print "%d time points for %d variables written to file %s" % (sol.ntimes, len(sol), fn)
        os.chdir(cwd)
    
    def orderByNames(self, varnames):
        for sol in self.solutions:
            sol.orderByNames(varnames)

    def orderByModelVars(self, amodel):
        vnames = [x for x in amodel().varnames]
        self.orderByNames(vnames)

def readTCs(source, filedir = None, intvarsorder = None, names = None, verbose = False):
    tcs = Solutions()
    tcsnames = None
    if isinstance(source, model.Model):
        #retrieve info from model declaration
        stcs = source['timecourses']
        tcs.filenames = stcs.filenames
        tcsnames = stcs.defaultnames
    else:
        tcs.filenames = source
    if names is None:
        if tcsnames is not None:
            names = tcsnames
    nread = tcs.loadTimeCourses(filedir, names=names, verbose=verbose)
    return tcs

TimeCourses = Solutions
Solution = SolutionTimeCourse

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
    demodata_noheader = """
#this is demo data without a header
#t x y z
0       1 0         0
0.1                  0.1

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.5  - 0.5 - -
0.6 0.6 0.8 0.9

"""

    
    aTC   = StringIO.StringIO(demodata)
    aTCnh = StringIO.StringIO(demodata_noheader)

    sol = SolutionTimeCourse()
    sol.load_from(aTC)   
    print '\n- using load_from() ----------------'
    print '\nnames:'
    print sol.names
    print '\nt'
    print sol.t
    print '\ndata'
    print sol.data
    print
    
    sol.load_from_str(demodata)
    sol.orderByNames("z y".split())
    print '\n- using load_from() with name order z y'
    print '\nnames:'
    print sol.names
    print '\ndata'
    print sol.data
    print

    sol.load_from_str(demodata)
    sol.orderByNames("z".split())
    print '\n- using load_from() with name order z'
    print '\nnames:'
    print sol.names
    print '\ndata'
    print sol.data
    print


    try:
        sol.load_from_str(demodata)
        print '\n- using load_from() with name order x bof z'
        sol.orderByNames("x bof z".split())
        print '\nnames:'
        print sol.names
        print '\ndata'
        print sol.data
        print
    except StimatorTCError, msg:
        print msg
        print

    sol.load_from_str(demodata)
    print '\n- using load_from() ----------------'
    print '\nnames:'
    print sol.names
    print '\nt'
    print sol.t
    print '\ndata'
    print sol.data
    print
    print '\n- now dumping, using save_to_str() ----------------'
    stc = sol.save_to_str()
    print stc
    print '-----------------------------------------------------'


    print '===Reading data without a header========================='
    aTCnh.seek(0)
    sol.load_from(aTCnh)   
    print '\n- using load_from(), names not provided'
    print '\nnames:'
    print sol.names
    print '\nt'
    print sol.t
    print '\ndata'
    print sol.data
    print
    aTCnh.seek(0)
    sol.load_from(aTCnh, names = ['v1','v2','v3', 'v4', 'v5'])   
    print '\n- using load_from() with names v1, v2 ,v3, v4, v5'
    print '\nnames:'
    print sol.names
    print '\nt'
    print sol.t
    print '\ndata'
    print sol.data
    print
    aTCnh.seek(0)
    sol.load_from(aTCnh, names = ['v1','v2'])   
    print '\n- using load_from() with names v1, v2'
    print '\nnames:'
    print sol.names
    print '\nt'
    print sol.t
    print '\ndata'
    print sol.data
    print
    
    #~ aTC.seek(0)
    #~ sol.load_from(aTC, atindexes=(0,3,1,2))   
    #~ print '\n- using load_from() atindexes (0,3,1,2)'
    #~ print '\nnames:'
    #~ print sol.names
    #~ print '\nt'
    #~ print sol.t
    #~ print '\ndata'
    #~ print sol.data
    #~ print
    
    print '==Using SolutionTimeCourse interface ===================='
    aTC.seek(0)
    sol.load_from(aTC)   
    print 'retrieving components...'
    try:
        print '\nnames:'
        print sol.names
        print '\nt'
        print sol.t
        print '\ndata'
        print sol.data
        print
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
        print 'for i in range(len(sol)): print sol[i]'
        for i in range(len(sol)):
            print sol[i]
        print 'for i in sol: print i'
        for i in sol:
            print i
        
        print "sol['k']"
        print sol['k']
    except ValueError, msg:
        print msg
    print
    print '\n- testing write_to() ----------------'
    sol.write_to('examples/exp.txt')
    print '\n- reading back from file ------------'
    sol.load_from('examples/exp.txt')
    print '\nnames:'
    print sol.names
    print '\nt'
    print sol.t
    print '\ndata'
    print sol.data
    print
    
    
    sol.load_from('examples/TSH2b.txt')
    print '\n- using load_from() ----------------'
    print '\nnames:'
    print sol.names
    print '\nnumber of times'
    print sol.ntimes
    print '\nshape'
    print sol.shape
    print '\nstate at 0.0:'
    print sol.state_at(0)
    print '\nlast time point:'
    print sol.last
    print
    
    sol2 = sol.clone()
    del(sol)
    print '\n- using a cloned solution ----------'
    print '\nnames:'
    print sol2.names
    print '\nnumber of times'
    print sol2.ntimes
    print '\nshape'
    print sol2.shape
    print '\nstate at 0.0:'
    print sol2.state_at(0)
    print '\nlast time point:'
    print sol2.last
    print
    
    sol = sol2.copy()
    del(sol2)
    print '\n- a cloned with copy() solution -----'
    print '\nnames:'
    print sol.names
    print '\nnumber of times'
    print sol.ntimes
    print '\nshape'
    print sol.shape
    print '\nstate at 0.0:'
    print sol.state_at(0)
    print '\nlast time point:'
    print sol.last
    print

    sol2 = sol.copy('HTA')
    del(sol)
    print "\n- a cloned with copy('HTA') solution --"
    print '\nnames:'
    print sol2.names
    print '\nnumber of times'
    print sol2.ntimes
    print '\nshape'
    print sol2.shape
    print '\nstate at 0.0:'
    print sol2.state_at(0)
    print '\nlast time point:'
    print sol2.last
    print

    print "-Reading tcs, using readTCs() -----------"
    tcs = readTCs(['TSH2b.txt', 'TSH2a.txt'], 'examples', verbose=True)
    for i, tc in enumerate(tcs):
        print tc.shape
        print tc.names
        print tc.state_at(0.0)
        print tc.last
        print tc.filename
        print tc.shortname
        print
    
    print "Providing default names HTA SDLTSH ------------------------"
    tcs = readTCs(['TSH2b.txt', 'TSH2a.txt'], 'examples', names = "SDLTSH HTA".split(), verbose=True)
    for i, tc in enumerate(tcs):
        print tc.shape
        print tc.names
        print tc.state_at(0.0)
        print tc.last
        print tc.shortname
        print
    
    print "After changing order to HTA SDLTSH ------------------------"
    
    tcs.orderByNames('HTA SDLTSH'.split())
    for i, tc in enumerate(tcs):
        print tc.shape
        print tc.names
        print tc.state_at(0.0)
        print tc.data[:,0]
        print tc.last
        print tc.shortname
        print
    
    print "saving to different files"
    tcs.saveTimeCoursesTo(['TSH2b_2.txt', 'TSH2a_2.txt'], 'examples', verbose=True)
    
    
    
    m = modelparser.read_model("""
    v1:        -> SDLTSH, rate = 1 ..
    v2: SDLTSH -> HTA,    rate = 2 ..
    timecourse TSH2b.txt
    timecourse TSH2a.txt
    variables SDLTSH HTA
    """)
    #print m
    
    print
    print
    print "After changing order according to model variables ------"
    
    tcs.orderByModelVars(m)
    for i, tc in enumerate(tcs):
        print tc.shape
        print tc.names
        print tc.state_at(0.0)
        print tc.data[:,0]
        print tc.last
        print tc.shortname
        print

    print "- Reading tcs using info declared in a model def -"
    tcs = readTCs(m, 'examples', verbose=True)
    for i, tc in enumerate(tcs):
        print tc.shape
        print tc.names
        print tc.state_at(0.0)
        print tc.last
        print tc.shortname
        print

    m = modelparser.read_model("""
    v1:        -> SDLTSH, rate = 1 ..
    v2: SDLTSH -> HTA,    rate = 2 ..
    timecourse ../stimator/examples/TSH2b.txt
    timecourse ../stimator/examples/TSH2a.txt
    variables SDLTSH HTA
    """)

    print "- Reading tcs using info declared in a model def -"
    print "(relative paths declared)"
    tcs = readTCs(m, verbose=True)
    for i, tc in enumerate(tcs):
        print tc.shape
        print tc.names
        print tc.state_at(0.0)
        print tc.last
        print tc.shortname
        print
