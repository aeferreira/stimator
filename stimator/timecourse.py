from __future__ import print_function, absolute_import, division
import os.path
from StringIO import StringIO
#from io import StringIO
import re
import numpy as np
import stimator.plots as plots
from stimator.utils import _is_string

FRAC_PATTERN = r"[-]?\d*[.]?\d+"
REAL_PATTERN = FRAC_PATTERN + r"(e[-]?\d+)?"
ID_RE = re.compile(r"[_a-z]\w*", re.IGNORECASE)
REAL_RE = re.compile(REAL_PATTERN, re.IGNORECASE)

class StimatorTCError(Exception):

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

# ----------------------------------------------------------------------------
#         THE BASIC TIMECOURSE CLASS
# ----------------------------------------------------------------------------

class SolutionTimeCourse(object):
    """Holds a timecourse created by ODE solvers"""

    def __init__(self, t=None, data=None, names=None, title="", dense=False):
        if t is None:
            t = np.array([])
        self.t = t          # values of time points
        if data is None:
            data = np.array([])
        self.data = data    # table of points: series in rows, times in cols
        if names is None:
            names = []
        self.names = names  # names of the series
        self.title = title  # a title for the solution
        self.dense = dense

        # for solutions read from a file
        self.filename = ""
        self.shortname = ""

    def __len__(self):
        """Retrieves the number of vars in this solution,
        NOT the len(timepoints)."""
        return self.data.shape[0]

    def __bool__(self):
        return len(self.t) > 0
    
    __nonzero__ = __bool__

    def __getNumberOfTimes(self):
        """Retrieves the number of time points"""
        return self.data.shape[1]
    ntimes = property(__getNumberOfTimes)

    def __getShape(self):
        return self.data.shape
    shape = property(__getShape)

    def __getitem__(self, key):
        """retrieves a series by name or index"""
        if _is_string(key):
            try:
                i = self.names.index(key)
            except ValueError:
                raise ValueError("No data for '%s' in timecourse" % str(key))
            return self.data.__getitem__(i)
        return self.data.__getitem__(key)

    def state_at(self, t):
        """Retrieves a dict with values at a time point.

           Interpolation may be necessary."""
        if t > self.t[-1] or t < self.t[0]:
            raise ValueError("No data for time '%s' in timecourse" % str(t))
        # Interpolate:
        ileft = self.t.searchsorted(t, side='left')
        iright = self.t.searchsorted(t, side='right')
        if iright == ileft:
            ileft -= 1
            tl = self.t[ileft]
            tr = self.t[iright]
            yl = self.data[:, ileft]
            yr = self.data[:, iright]
            m = (yr - yl) / (tr - tl)
            y = yl + m * (t - tl)
        else:
            y = self.data[:, ileft]
        return dict([(x, value) for (x, value) in zip(self.names, y)])

    def i_time(self, t):
        """Retrieves the closest index for time t."""
        if t > self.t[-1] or t < self.t[0]:
            raise ValueError("No data for time '%s' in timecourse" % str(t))
        # Find closest:
        ileft = self.t.searchsorted(t, side='left')
        iright = self.t.searchsorted(t, side='right')
        if iright == ileft:
            ileft -= 1
            tl = self.t[ileft]
            tr = self.t[iright]
            if (t - tl) <= (tr - t):
                return ileft
            else:
                return iright
        else:
            return ileft

    def __getLastState(self):
        """Retrieves state_at last timepoint"""
        return self.state_at(self.t[-1])
    last = property(__getLastState)  # use as 'sol.last'

    def __getInitState(self):
        """Retrieves state_at first timepoint"""
        return self.state_at(self.t[0])
    init = property(__getInitState)  # use as 'sol.init'

    def apply_transf(self, f, newnames=None, new_title=None):
        """Apply a transformation to time series in place.

           f is the transformation function, with signature
           f(variables,t). variables is an array, list or tuple, t is a scalar.
           newnames is a list of names of the transformed variables.
           results are kept 'in place': data is substituted."""

        def newf(newdata, f):
            return f(newdata[1:], newdata[0])
        trf = np.apply_along_axis(newf, 0, np.vstack((self.t, self.data)), f)
        if newnames is not None:
            self.names = newnames
        if new_title is not None:
            self.title = new_title
        self.data = trf
        return self

    def transform(self, f, newnames=None, new_title=None):
        """Apply a transformation to time series.

           f is the transformation function, with signature
           f(variables,t). variables is an array, list or tuple, t is a scalar.
           newnames is a list of names of the transformed variables."""

        return self.clone().apply_transf(f, newnames, new_title)

    def read_str(self, s, names=None):
        aTC = StringIO(s)
        aTC.seek(0)
        result = self.read_from(aTC, names)
        aTC.close()
        return result

    def read_from(self, filename, names=None):
        """Reads a time course from file.

        Fills self.names from a header with variable names
        (possibly absent in file). Fills a 2D numpy array with
        whitespace separated data.
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
            f = open(filename, "rU")  # could be a name,instead of an open file
            isname = True

        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue    # empty lines are skipped
            if line.startswith('#'):
                continue    # comment lines are skipped
            items = line.split()

            if ID_RE.match(items[0]):
                if not headerFound and not t0found:
                    header = [item for item in items if ID_RE.match(item)]
                    headerFound = True
                else:
                    continue
            elif not REAL_RE.match(items[0]):
                continue
            else:
                if not t0found:
                    nvars = len(items)
                    t0found = True
                temprow = [np.nan] * nvars
                for (i, num) in enumerate(items):
                    if REAL_RE.match(num):
                        temprow[i] = float(num)
                rows.append(temprow)
        if isname:
            f.close()

        # create default names "t, x1, x2, x3,..." or use names if provided
        if len(header) == 0:
            header = ['t']
            for i in range(1, nvars):
                header.append('x%d' % i)
            if names is not None:
                smallindx = min(len(header) - 1, len(names))
                for i in range(smallindx):
                    header[i + 1] = names[i]
        data = np.array(rows)
        self.names = header[1:]
        self.t = data[:, 0].T
        self.data = data[:, 1:].T
        return self

    def _write_to_str(self):
        aTC = StringIO()
        aTC.seek(0)
        self.write_to(aTC)
        return aTC.getvalue()

    def __str__(self):
        return self._write_to_str()

    def write_to(self, filename):
        """Writes a time course to a file or file-like object.
        """

        if hasattr(filename, 'read'):
            f = filename
            isname = False
        else:
            f = open(filename, "w")  # could be a name,instead of an open file
            isname = True

        f.write("%s %s\n" % ('t', " ".join(self.names)))
        npoints = len(self.t)
        for i in range(npoints):
            row = [self.t[i]]
            row.extend(self.data[:, i])
            row = " ".join([str(j) for j in row])
            f.write("%s\n" % row)
        if isname:
            f.close()

    def clone(self, new_title=None):
        """Clones the entire solution."""
        tc = SolutionTimeCourse(self.t.copy(),
                                self.data.copy(),
                                self.names[:],
                                self.title, self.dense)
        if new_title is not None:
            tc.title = new_title
        tc.filename = self.filename
        tc.shortname = self.shortname
        return tc

    def copy(self, names=[], newtitle=None):
        """Constructs new solution, restricted to the variables in 'names'."""
        if not (isinstance(names, list) or isinstance(names, tuple)):
            names = names.strip()
            names = names.split()
        t = self.t.copy()
        if names == []:
            names = self.names
        nameindexes = []
        for name in names:
            if name not in self.names:
                raise ValueError("No data for '%s' in timecourse" % name)
            nameindexes.append(self.names.index(name))
        data = self.data[nameindexes, :].copy()
        if newtitle is not None:
            title = newtitle
        else:
            title = self.title
        tc = SolutionTimeCourse(t, data, names[:], title, dense=self.dense)
        tc.filename = self.filename
        tc.shortname = self.shortname
        return tc

    def orderByNames(self, varnames):
        oldindexes = list(range(len(self)))
        newindexes = []
        for vname in varnames:
            if vname in self.names:
                indx = self.names.index(vname)
                newindexes.append(indx)
                oldindexes.remove(indx)
        newindexes.extend(oldindexes)
        self.names = [self.names[i] for i in newindexes]
        self.data = self.data[np.array(newindexes, dtype=int)]

    def plot(self, axes=None, **kwargs):
        if axes is not None:
            axis_set = [axes]
        else:
            axis_set = None
        ss = Solutions([self])
        ss.plot(axis_set=axis_set, **kwargs)

# ----------------------------------------------------------------------------
#         A CONTAINER FOR SEVERAL TIMECOURSES
# ----------------------------------------------------------------------------


class Solutions(object):
    """Holds a colection of objects of class SolutionTimeCourse"""

    def __init__(self, aList=None, title=""):
        self.title = title
        self.solutions = []
        self.shortnames = []
        self.filenames = []
        self.basedir = None
        self.defaultnames = None  # list of names to use if headers are missing
        # aList argument must be an iterable
        # TODO: throw Exception if it isn't
        if aList is not None:
            for s in aList:
                self.append(s)

    def __str__(self):
        output = (str(s) for s in self.solutions)
        return '\n'.join(output)

    def __getitem__(self, key):
        """retrieves a series by index"""
        return self.solutions.__getitem__(key)

    def __len__(self):
        return len(self.solutions)

    def __bool__(self):
        return len(self.solutions) > 0
    
    __nonzero__=__bool__

    def __iadd__(self, other):
        if isinstance(other, Solutions):
            self.solutions.extend(other.solutions)
        elif isinstance(other, list) or isinstance(other, tuple):
            for s in other:
                if not isinstance(s, SolutionTimeCourse):
                    raise TypeError("Must add a solution or a set of them")
            self.solutions.extend(list(other))
        elif isinstance(other, SolutionTimeCourse):
            self.solutions.append(other)
        else:
            raise TypeError("Must add a solution or a set of solutions")
        return self

    def __iter__(self):
        return iter(self.solutions)

    def append(self, other):
        return self.__iadd__(other)

    def loadTimeCourses(self, filedir=None, names=None, verbose=False):
        if len(self.filenames) == 0:
            msg = "No time courses to load!"
            msg += "Please indicate time courses with 'timecourse <filename>'"
            raise StimatorTCError(msg)

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
            print ("-- reading time courses -------------------------------")
        for filename in pathlist:
            if not os.path.exists(filename) or not os.path.isfile(filename):
                os.chdir(cwd)
                raise StimatorTCError("File \n%s\ndoes not exist" % filename)
            sol = SolutionTimeCourse()
            sol.read_from(filename, names=names)
            if sol.shape == (0, 0):
                os.chdir(cwd)
                error_msg = "File\n%s\ndoes not contain valid data" % filename
                raise StimatorTCError(error_msg)
            else:
                if verbose:
                    print("file %s:" % (filename))
                    print("%d time points, %d variables" % (sol.ntimes,
                                                            len(sol)))
                self.append(sol)
                nTCsOK += 1
        self.shortnames = [os.path.split(filename)[1] for filename in pathlist]
        for i, sol in enumerate(self.solutions):
            sol.title = self.shortnames[i]
            sol.shortname = self.shortnames[i]
            sol.filename = self.filenames[i]
        os.chdir(cwd)
        return nTCsOK

    def saveTimeCoursesTo(self, filenames, filedir=None, verbose=False):
        if len(self) == 0:
            print ("No time courses to save!")
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
            print ("-------------------------------------------------------")
        for fn, sol in zip(pathlist, self.solutions):
            sol.write_to(fn)
            if verbose:
                print ("%d time points for %d variables written to file %s" % (sol.ntimes, len(sol), fn))
        os.chdir(cwd)

    def orderByNames(self, varnames):
        for sol in self.solutions:
            sol.orderByNames(varnames)
        return self

    def orderByModelVars(self, amodel):
        vnames = [x for x in amodel.varnames]
        self.orderByNames(vnames)
        return self

    def plot(self, **kwargs):
        return plots.plotTCs(self, **kwargs)
    

def readTCs(source,
            filedir=None,
            intvarsorder=None,
            names=None,
            verbose=False):
    tcs = Solutions()
    tcsnames = None
    if hasattr(source, 'metadata'):
        # retrieve info from model declaration
        stcs = source.metadata['timecourses']
        tcs.filenames = stcs.filenames
        tcsnames = stcs.defaultnames
    elif _is_string(source):
        tcs.filenames = [source]
    else:
        tcs.filenames = source
    if names is None:
        if tcsnames is not None:
            names = tcsnames
    tcs.loadTimeCourses(filedir, names=names, verbose=verbose)
    return tcs

read_tc = readTCs
TimeCourses = Solutions
Solution = SolutionTimeCourse

# ----------------------------------------------------------------------------
#         Time course divergence metrics
# ----------------------------------------------------------------------------

def extendedKLdivergence(modelTCs, deltaT, indexes):
    result = []
    for (i, j) in indexes:
        m = modelTCs[i].data
        n = modelTCs[j].data
        m = np.where(m <= 0.0, np.NaN, m)
        n = np.where(n <= 0.0, np.NaN, n)
        dif = -deltaT * np.nansum(np.float64(m * (np.log(m / n) + n / m - 1.0)))
        result.append(dif)
    return result


def KLdivergence(modelTCs, deltaT, indexes):
    result = []
    for (i, j) in indexes:
        m = modelTCs[i].data
        n = modelTCs[j].data
        m = np.where(m <= 0.0, np.NaN, m)
        n = np.where(n <= 0.0, np.NaN, n)
        dif = -deltaT * np.nansum(np.float64(m * np.log(m / n)))
        result.append(dif)
    return result


def L2_midpoint_weights(modelTCs, deltaT, indexes):
    """L2-norm for time courses, weighted by midpoints"""

    result = []
    for i in range(len(modelTCs) - 1):
        for j in range(i + 1, len(modelTCs)):
            numResult = 0.0
            for tc1, tc2 in zip(modelTCs[i], modelTCs[j]):
                tempTC = np.float64((((tc1 - tc2)**2) / (((tc1 + tc2)/2.0)**2)) * deltaT)
                numResult -= np.nansum(tempTC)
            result.append(numResult)
    return result


def L2(modelTCs, deltaT, indexes):
    """L2-norm for time courses"""
    result = []
    for i in range(len(modelTCs) - 1):
        for j in range(i + 1, len(modelTCs)):
            numResult = 0.0
            for tc1, tc2 in zip(modelTCs[i], modelTCs[j]):
                tempTC = np.float64(((tc1 - tc2) ** 2)) * deltaT
                numResult -= np.nansum(tempTC)
            result.append(numResult)
    return result


def _transform2array(vect):
    if isinstance(vect, float) or isinstance(vect, int):
        res = np.array((vect), dtype=float)
    elif isinstance(vect, list) or isinstance(vect, tuple):
        res = np.diag(np.array(vect, dtype=float))
    else:
        res = vect  # is already an array (must be 2D)
    return res


def constError_func(vect):
    res = _transform2array(vect)

    def CE(x):
        return res
    return CE


def propError_func(vect):
    res = _transform2array(vect)

    def CE(x):
        return res * x
    return CE


def getFullTCvarIndexes(model, tcs):
    # mask series with NaN values.
    allmodelvarindexes, alltcvarindexes = [], []
    for data in tcs:
        nt = data.ntimes
        varindexes = []
        modelvarindexes = []

        for ivar in range(len(data.data)):
            # count NaN
            yexp = data[ivar]
            nnan = len(yexp[np.isnan(yexp)])
            if nnan >= nt - 1:
                continue
            varindexes.append(ivar)
            vname = data.names[ivar]
            indx = model.varnames.index(vname)
            modelvarindexes.append(indx)
        alltcvarindexes.append(np.array(varindexes, int))
        allmodelvarindexes.append(np.array(modelvarindexes, int))
    return allmodelvarindexes, alltcvarindexes


def getCommonFullVars(tcs):
    """Returns a list of names of variables that
    have full data in all timecourses."""

    common_names = []
    for itc, tc in enumerate(tcs):
        nt = tc.ntimes
        tcnames = tc.names
        for i, line in enumerate(tc.data):
            # count NaN
            yexp = line
            xname = tcnames[i]
            nnan = len(yexp[np.isnan(yexp)])
            if nnan >= nt - 1:
                if xname in common_names:
                    common_names.remove(xname)
            else:
                if itc == 0:
                    common_names.append(xname)
    return common_names


def getRangeVars(tcs, varnames):
    ranges = [0.0 for i in range(len(varnames))]
    for ix, x in enumerate(varnames):
        for tc in tcs:
            yexp = tc[x]
            tpe = (max(yexp) - min(yexp))
            ranges[ix] = max(ranges[ix], tpe)
    return ranges


def getCriteriumFunction(weights, model, tc):
    """Returns a function to compute the objective function
    (for each timecourse).

    The function has signature
    criterium(Y,i)
    Y is the predicted timecourse, for a given set of parameters.
    i is the index of the timecourse.
    The function returns a float.

    tc is a Solutions object holding ('experimental') timecourse data,
    each timecourse has shape (nvars, ntimes).

    weights can be:

    None         : no weighting (simple LSquares, S = sum((Ypred-Yexp)**2))
    all others are weighted LSquares, S = (Ypred-Yexp).T * W * (Ypred-Yexp)
    'demo'       : demo weighting  W = 1/j with j = 1,...,nvars
    """

    allmodelvarindexes, alltcvarindexes = getFullTCvarIndexes(model, tc)

    if weights is None:
        def criterium(Y, i):
            d = (Y.T[allmodelvarindexes[i]] - tc[i].data[alltcvarindexes[i]])
            return np.sum(d * d)
        return criterium

    if weights == 'demo':
        W = []
        for i in range(len(tc)):
            W.append(np.array([1.0 / (1 + j) for j in range(alltcvarindexes[i])]))

        # print W
        def criterium(Y, i):
            d = (Y.T[allmodelvarindexes[i]] - tc.data[i][alltcvarindexes[i]])
            return np.sum(d * W[i] * d)
        return criterium

    # TODO: weights not implemented
    return None

# ----------------------------------------------------------------------------
#         TESTING CODE
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    from modelparser import read_model
    print ('\n===Parsing in-code timecourse ========================')

    demodata = """
#this is demo data with a header
t x y z
0       0.95 0         0
0.1                  0.1

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.3 0.4 0.5 0.55
0.4 0.5 0.6 0.7
0.5 0.6 0.8 0.9
0.55 0.7 0.85 0.95
0.6  - 0.5 - -

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

    demodata2 = """
#this is demo data with a header
t x y z
0       0.95 0         0
0.1                  0.09

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.3 0.45 0.55 0.58
0.4 0.5 0.65 0.75
0.5 0.65 0.85 0.98
0.55 0.7 0.9 1.45
0.6  - 0.4 - -
"""

    aTC = StringIO(demodata)

    sol = SolutionTimeCourse().read_from(aTC)
    print ('\n!! using read_from() ----------------')
    print ('\nnames:')
    print (sol.names)
    print ('\nt')
    print (sol.t)
    print ('\ndata')
    print (sol.data, '\n')
    print('\n!! printing a solution ----------------')
    print (sol)
    print('----------------------------------------')

    sol = SolutionTimeCourse().read_str(demodata)
    sol.orderByNames("z y".split())
    print ('\n!! using read_str() with name order z y')
    print ('\nnames:')
    print (sol.names)
    print ('\ndata')
    print (sol.data, '\n')
    
    sol.read_str(demodata)
    sol.orderByNames("z".split())
    print ('\n!! using read_str() with name order z')
    print ('\nnames:')
    print (sol.names)
    print ('\ndata')
    print (sol.data, '\n')

    try:
        sol.read_str(demodata)
        print ('\n!! using read_str() with name order x bof z')
        sol.orderByNames("x bof z".split())
        print ('\nnames:')
        print (sol.names)
        print ('\ndata')
        print (sol.data, '\n')
    except StimatorTCError as msg:
        print (msg, '\n')

    sol.read_str(demodata)
    print ('\n!! using read() ----------------')
    print ('\nnames:')
    print (sol.names)
    print ('\nt')
    print (sol.t)
    print ('\ndata')
    print (sol.data, '\n')

    print ('\n!! displaying with print() ----------------')
    # stc = str(sol)
    print (sol)
    print ('-----------------------------------------------------')

    print ('===Reading data without a header=========================')
    sol.read_str(demodata_noheader)
    print ('\n!! using read_from(), names not provided')
    print ('\nnames:')
    print (sol.names)
    print ('\nt')
    print (sol.t)
    print ('\ndata')
    print (sol.data, '\n')

    names = ['v1', 'v2', 'v3', 'v4', 'v5']
    sol.read_str(demodata_noheader, names=names)
    print ('\n!! using read_from() with names v1, v2 ,v3, v4, v5')
    print ('\nnames:')
    print (sol.names)
    print ('\nt')
    print (sol.t)
    print ('\ndata')
    print (sol.data, '\n')

    sol.read_str(demodata_noheader, names=['v1', 'v2'])
    print ('\n!! using read_from() with names v1, v2')
    print ('\nnames:')
    print (sol.names)
    print ('\nt')
    print (sol.t)
    print ('\ndata')
    print (sol.data, '\n')

    print ('==Using SolutionTimeCourse interface ====================')
    sol.read_str(demodata)
    print ('retrieving components...')
    try:
        print ('\nnames:')
        print (sol.names)
        print ('\nt')
        print (sol.t)
        print ('\ndata')
        print (sol.data, '\n')

        print ('len(sol)')
        print (len(sol))
        print ('sol.ntimes')
        print (sol.ntimes)
        print ('sol[0] (first var, "x")')
        print (sol[0])
        print ('sol.t')
        print (sol.t)
        print ("sol['x']")
        print (sol['x'])
        print ("sol.names")
        print (sol.names)
        print ('Last time point, sol[:,-1] returns array')
        print (sol[:, -1])
        print (type(sol[:, -1]))
        print ('The following return dictionaries:')
        print ('sol.state_at(0.2)')
        print (sol.state_at(0.2))
        print (type(sol.state_at(0.2)))
        print ('sol.state_at(0.55)')
        print (sol.state_at(0.55))
        print ('sol.state_at(0.0)')
        print (sol.state_at(0.0))
        print ('sol.state_at(0.6)')
        print (sol.state_at(0.6))
        print ('sol.last (Last time point the easy way)')
        print (sol.last)
        print (type(sol.last))
        print ('sol.last["x"]')
        print (sol.last['x'])
        print ('sol.init (Initial point the easy way)')
        print (sol.init)
        print (type(sol.init))
        print ('sol.init["x"]')
        print (sol.init['x'])
        print ('for i in range(len(sol)): print sol[i]')
        for i in range(len(sol)):
            print (sol[i])
        print ('for i in sol: print i')
        for i in sol:
            print (i)

        print ("sol['k']")
        print (sol['k'])
    except ValueError as msg:
        print (msg)
    print ('\n')
    print ('\n!! testing write_to() ----------------')
    sol.write_to('examples/exp.txt')
    print ('\n!! reading back from file ------------')
    sol.read_from('examples/exp.txt')
    print ('\nnames:')
    print (sol.names)
    print ('\nt')
    print (sol.t)
    print ('\ndata')
    print (sol.data, '\n')

    print ('== Using Solutions interface ====================')

    sols = Solutions(title='all time courses')

    print ('\n!! testing boolean context on Solutions ----------------')
    
    if sols:
        print ('Sols has solutions')
    else:
        print ('Sols is empty')

    s = SolutionTimeCourse(title='the first time course').read_str(demodata)
    sols += s
    s = SolutionTimeCourse(title='the second time course').read_str(demodata2)
    sols += s

    if sols:
        print ('Sols has solutions')
    else:
        print ('Sols is empty')
    
    print ('\n!! printing Solutions ----------------')
    
    print(sols)

    print ('\n!! testing transformations ----------------')

    sols = Solutions(title='all time courses')

    s = SolutionTimeCourse(title='original time course').read_str(demodata2)
    sols += s
    print ('--- before transformation')
    print ('- names')
    print (s.names)
    print ('- data')
    print (s.data)

    def average(x, t):
        # print ('applying transformation')
        return np.array([t/2.0, (x[0]+x[-1])/2.0])

    s = s.transform(average,
                    newnames=['t/2', 'mid point'],
                    new_title='after transformation')
    print ('--- after transformation')
    print ('- names')
    print (s.names)
    print ('- data')
    print (s.data)
    sols += s

    sol.read_from('examples/timecourses/TSH2b.txt')
    print ('\n!! using read() ----------------')
    print ('\nnames:')
    print (sol.names)
    print ('\nnumber of times')
    print (sol.ntimes)
    print ('\nshape')
    print (sol.shape)
    print ('\nstate at 0.0:')
    print (sol.state_at(0))
    print ('\nlast time point:')
    print (sol.last, '\n')

    sol2 = sol.clone()
    del sol
    print ('\n!! using a cloned solution ----------')
    print ('\nnames:')
    print (sol2.names)
    print ('\nnumber of times')
    print (sol2.ntimes)
    print ('\nshape')
    print (sol2.shape)
    print ('\nstate at 0.0:')
    print (sol2.state_at(0))
    print ('\nlast time point:')
    print (sol2.last, '\n')

    sol = sol2.copy()
    del sol2
    print ('\n!! a cloned with copy() solution -----')
    print ('\nnames:')
    print (sol.names)
    print ('\nnumber of times')
    print (sol.ntimes)
    print ('\nshape')
    print (sol.shape)
    print ('\nstate at 0.0:')
    print (sol.state_at(0))
    print ('\nlast time point:')
    print (sol.last, '\n')

    sol2 = sol.copy('HTA')
    del sol
    print ("\n!! a cloned with copy('HTA') solution --")
    print ('\nnames:')
    print (sol2.names)
    print ('\nnumber of times')
    print (sol2.ntimes)
    print ('\nshape')
    print (sol2.shape)
    print ('\nstate at 0.0:')
    print (sol2.state_at(0))
    print ('\nlast time point:')
    print (sol2.last, '\n')

    print ("-Reading tcs, using readTCs() -----------")
    tcs = readTCs(['TSH2b.txt', 'TSH2a.txt'], 'examples/timecourses',
                  verbose=True)
    print ('\nResulting type:', type(tcs))
    print ('\nElements:')
    for i, tc in enumerate(tcs):
        print (i, '---->>>')
        print ('type:', type(tc))
        print ('shape', tc.shape)
        print ('names:', tc.names)
        print ('state at 0.0', tc.state_at(0.0))
        print ('last', tc.last)
        print ('filename:', tc.filename)
        print ('shortname:', tc.shortname, '\n')

    print ("-Providing default names HTA SDLTSH ------------------------")
    tcs = readTCs(['TSH2b.txt', 'TSH2a.txt'],
                  'examples/timecourses',
                  names="SDLTSH HTA".split(),
                  verbose=True)
    print ('\nResulting type:', type(tcs))
    print ('\nElements:')
    for i, tc in enumerate(tcs):
        print (i, '---->>>')
        print ('type:', type(tc))
        print ('shape', tc.shape)
        print ('names:', tc.names)
        print ('state at 0.0', tc.state_at(0.0))
        print ('last', tc.last)
        print ('filename:', tc.filename)
        print ('shortname:', tc.shortname, '\n')

    print ("After changing order to HTA SDLTSH ------------------------")

    tcs.orderByNames('HTA SDLTSH'.split())
    for i, tc in enumerate(tcs):
        print (i, '---->>>')
        print ('type:', type(tc))
        print ('shape', tc.shape)
        print ('names:', tc.names)
        print ('data[:,0] ->', tc.data[:, 0])
        print ('last', tc.last)
        print ('shortname:', tc.shortname, '\n')

    print ("saving to different files")
    tcs.saveTimeCoursesTo(['TSH2b_2.txt', 'TSH2a_2.txt'],
                          'examples',
                          verbose=True)

    m = read_model("""
    v1:        -> SDLTSH, rate = 1 ..
    v2: SDLTSH -> HTA,    rate = 2 ..
    timecourse TSH2b.txt
    timecourse TSH2a.txt
    variables SDLTSH HTA
    """)
    # print m

    print ('\n\n')
    print ("After changing order according to model variables ------")

    tcs.orderByModelVars(m)
    for i, tc in enumerate(tcs):
        print (i, '---->>>')
        print ('type:', type(tc))
        print ('shape', tc.shape)
        print ('names:', tc.names)
        print ('data[:,0] ->', tc.data[:, 0])
        print ('last', tc.last)
        print ('shortname:', tc.shortname, '\n')

    print ("!! Reading tcs using info declared in a model def -")
    tcs = readTCs(m, 'examples/timecourses', verbose=True)
    print ('\nResulting type:', type(tcs))
    print ('\nElements:')
    for i, tc in enumerate(tcs):
        print (i, '---->>>')
        print ('type:', type(tc))
        print ('shape', tc.shape)
        print ('names:', tc.names)
        print ('state at 0.0', tc.state_at(0.0))
        print ('last', tc.last)
        print ('filename:', tc.filename)
        print ('shortname:', tc.shortname, '\n')

    m = read_model("""
    v1:        -> SDLTSH, rate = 1 ..
    v2: SDLTSH -> HTA,    rate = 2 ..
    timecourse ../stimator/examples/timecourses/TSH2b.txt
    # timecourse ../stimator/examples/timecourses/TSH2a.txt
    variables SDLTSH HTA
    """)

    print ("!! Reading tcs using info declared in a model def -")
    print ("(relative paths declared)")
    tcs = readTCs(m, verbose=True)
    print ('\nResulting type:', type(tcs))
    print ('\nElements:')
    for i, tc in enumerate(tcs):
        print (i, '---->>>')
        print ('type:', type(tc))
        print ('shape', tc.shape)
        print ('names:', tc.names)
        print ('state at 0.0', tc.state_at(0.0))
        print ('last', tc.last)
        print ('filename:', tc.filename)
        print ('shortname:', tc.shortname, '\n')
