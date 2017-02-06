from __future__ import print_function, absolute_import, division
import os.path
import re
import numpy as np

from six import StringIO
from six.moves import getcwd

import stimator.plots as plots
from stimator.utils import _is_string, _is_number, _is_sequence

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

    def read_from(self, source, names=None):
        """Reads a time course from a path or file-like object.

        Fills self.names from a header with variable names
        (possibly absent in file). Fills a 2D numpy array with
        whitespace separated data.
        """

        header = []
        nvars = 0
        rows = []
        headerFound = False
        t0found = False

        ishandle = False
        try:
            f = open(source, "rU")
        except TypeError:
            ishandle = True
            f = source

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
        if not ishandle:
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

    def __str__(self):
        out = ["%s %s" % ('t', " ".join(self.names))]
        npoints = len(self.t)
        for i in range(npoints):
            row = [self.t[i]]
            row.extend(self.data[:, i])
            row = " ".join([str(j) for j in row])
            out.append(row)
        return '\n'.join(out)

    def write_to(self, destination):
        """Writes a time course to a path or file-like object.
        """

        ishandle = False
        try:
            f = open(destination, "w")
        except TypeError:
            ishandle = True
            f = destination
        
        f.write(str(self))
        if not ishandle:
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

    def order_by_names(self, varnames):
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

    def _load_tcs(self, filedir=None, names=None, verbose=False):
        if len(self.filenames) == 0:
            msg = "No time courses to load!"
            msg += "Please indicate time courses with 'timecourse <filename>'"
            raise StimatorTCError(msg)

        # check and load timecourses
        cwd = getcwd()
        if filedir is None:
            filedir = ''
        
        plist = [os.path.join(cwd, filedir, k) for k in self.filenames]
        plist = [os.path.abspath(p) for p in plist]

        self.data = []
        nTCsOK = 0
        if verbose:
            print ("-- reading time courses -------------------------------")
        for filename in plist:
            if not os.path.exists(filename) or not os.path.isfile(filename):
                raise StimatorTCError("File \n%s\ndoes not exist" % filename)
            sol = SolutionTimeCourse()
            sol.read_from(filename, names=names)
            if sol.shape == (0, 0):
                error_msg = "File\n%s\ndoes not contain valid data" % filename
                raise StimatorTCError(error_msg)
            else:
                if verbose:
                    print("file %s:" % (filename))
                    print("%d time points, %d variables" % (sol.ntimes,
                                                            len(sol)))
                self.append(sol)
                nTCsOK += 1
        self.shortnames = [os.path.split(filename)[1] for filename in plist]
        for i, sol in enumerate(self.solutions):
            sol.title = self.shortnames[i]
            sol.shortname = self.shortnames[i]
            sol.filename = self.filenames[i]
        return nTCsOK

    def write_to(self, filenames, filedir=None, verbose=False):
        if len(self) == 0:
            print ("No time courses to save!")
            return 0

        # check and load timecourses
        cwd = getcwd()
        if filedir is None:
            filedir = ''
        names = [os.path.join(cwd, filedir, k) for k in filenames]
        names = [os.path.abspath(p) for p in names]

        fstring = "{} time points for {} variables written to file {}".format

        if verbose:
            print ("-------------------------------------------------------")
        for name, sol in zip(names, self.solutions):
            sol.write_to(name)
            if verbose:
                print(fstring(sol.ntimes, len(sol), name))

    def order_by_names(self, varnames):
        for sol in self.solutions:
            sol.order_by_names(varnames)
        return self

    def order_by_modelvars(self, amodel):
        vnames = [x for x in amodel.varnames]
        self.order_by_names(vnames)
        return self
    
    def get_common_full_vars(self):
        """Return a list of names of variables that
        have full data in all timecourses."""

        for itc, tc in enumerate(self):
            names = []
            nt = tc.ntimes
            for name, yexp in zip(tc.names, tc.data):
                nnan = len(yexp[np.isnan(yexp)])
                if nnan < nt - 1:
                    names.append(name)
            if itc == 0:
                common_names = set(names)
            else:
                common_names.intersection_update(names)
        return list(common_names)

    def plot(self, **kwargs):
        return plots.plotTCs(self, **kwargs)
    

def read_tc(source,
            filedir=None,
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
    tcs._load_tcs(filedir, names=names, verbose=verbose)
    return tcs

# convenient or backwards compatible aliases
readTCs = read_tc
TimeCourses = Solutions
Solution = SolutionTimeCourse

# ----------------------------------------------------------------------------
#         Time course divergence metrics
# ----------------------------------------------------------------------------

def extendedKLdivergence(timecourses, delta_t, index_list):
    result = []
    for i, j in index_list:
        m = timecourses[i].data
        n = timecourses[j].data
        m = np.where(m <= 0.0, np.NaN, m)
        n = np.where(n <= 0.0, np.NaN, n)
        d = -delta_t * np.nansum(np.float64(m * (np.log(m / n) + n / m - 1.0)))
        result.append(d)
    return result


def KLdivergence(timecourses, delta_t, index_list):
    result = []
    for i, j in index_list:
        m = timecourses[i].data
        n = timecourses[j].data
        m = np.where(m <= 0.0, np.NaN, m)
        n = np.where(n <= 0.0, np.NaN, n)
        d = -delta_t * np.nansum(np.float64(m * np.log(m / n)))
        result.append(d)
    return result


def L2_midpoint_weights(timecourses, delta_t, indexes):
    """L2-norm for time courses, weighted by midpoints"""

    result = []
    for i in range(len(timecourses) - 1):
        for j in range(i + 1, len(timecourses)):
            numResult = 0.0
            for tc1, tc2 in zip(timecourses[i], timecourses[j]):
                tempTC = np.float64((((tc1 - tc2)**2) / (((tc1 + tc2)/2.0)**2)) * delta_t)
                numResult -= np.nansum(tempTC)
            result.append(numResult)
    return result


def L2(timecourses, delta_t, indexes):
    """L2-norm for time courses"""
    result = []
    for i in range(len(timecourses) - 1):
        for j in range(i + 1, len(timecourses)):
            numResult = 0.0
            for tc1, tc2 in zip(timecourses[i], timecourses[j]):
                tempTC = np.float64(((tc1 - tc2) ** 2)) * delta_t
                numResult -= np.nansum(tempTC)
            result.append(numResult)
    return result


def _transform2array(vect):
    """Given a float or sequence, transform into a diagonal array.
       A 2D array is left unchanged."""
    if _is_number(vect):
        res = np.array((vect), dtype=float)
    elif _is_sequence(vect):
        res = np.diag(np.array(vect, dtype=float))
    else:
        res = vect # must be already an array (must be 2D)
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


def get_fulltc_indexes(model, tcs):
    # mask series with NaN values.
    all_model_indexes, all_tc_indexes = [], []
    for tc in tcs:
        ntimes = tc.ntimes
        vars_i = []
        vars_i_model = []

        for i, yexp in enumerate(tc):
            # count NaN
            nnan = len(yexp[np.isnan(yexp)])
            if nnan >= ntimes - 1:
                continue
            vars_i.append(i)
            vars_i_model.append(model.varnames.index(tc.names[i]))
        
        all_tc_indexes.append(np.array(vars_i, int))
        all_model_indexes.append(np.array(vars_i_model, int))
    return all_model_indexes, all_tc_indexes


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

    vars_i_model, vars_i_tcs = get_fulltc_indexes(model, tc)

    if weights is None:
        def criterium(Y, i):
            d = (Y.T[vars_i_model[i]] - tc[i].data[vars_i_tcs[i]])
            return np.sum(d * d)
        return criterium

    if weights == 'demo':
        W = []
        for i in range(len(tc)):
            W.append(np.array([1.0 / (1 + j) for j in range(vars_i_tcs[i])]))

        # print W
        def criterium(Y, i):
            d = (Y.T[vars_i_model[i]] - tc.data[i][vars_i_tcs[i]])
            return np.sum(d * W[i] * d)
        return criterium

    # TODO: weights not implemented
    return None
