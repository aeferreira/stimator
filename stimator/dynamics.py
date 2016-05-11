from __future__ import print_function, division, absolute_import
import re
import math
#import itertools
from six.moves import zip as izip
from numpy import *
from scipy import integrate
from stimator.timecourse import SolutionTimeCourse, Solutions

from stimator.examples import models

identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)

def identifiersInExpr(_expr):
    iterator = identifier.finditer(_expr)
    return [_expr[m.span()[0]:m.span()[1]] for m in iterator]


def _is_string(a):
    return (isinstance(a, str) or
            isinstance(a, unicode))

def _is_sequence(arg):
    return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))


def _find_indexof_component(model, name):
    for i, p in enumerate(model.with_bounds):
        if p.name == name:
            return i, 'uncertain'
    for p in model.parameters:
        if p.name == name:
            return -1, 'parameters'
    i = model._Model__reactions.iget(name)
    if i is not None:
        return i, 'reactions'
    try:
        i = model._Model__variables.index(name)
        return i, 'variables'
    except:
        pass
    i = model._Model__transf.iget(name)
    if i is not None:
        return i, 'transf'
    i = model._Model__invars.iget(name)
    if i is not None:
        return i, 'invar'
    raise AttributeError('%s is not a component in this model' % name)



def init2array(model):
    """Transforms a state object into a numpy.array object.
       
       This is necessary for most numerical functions of numpy+scipy.
       Can accept the name of a state (must exist in Model) or state object.
       Values are returned in the order of model variables.
    """
    return array([model.get_init(var) for var in model.varnames])

def genStoichiometryMatrix(m):
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)

    vnames = m.varnames
    N = zeros((len(vnames), len(m.reactions)), dtype=float)
    for j, v in enumerate(m.reactions):
        for rORp, sign_one in [(v._reagents, -1.0),(v._products, 1.0)]:
            for var, coef in rORp:
                if var in vnames:
                    ivar = vnames.index(var)
                    N[ivar, j] = coef * sign_one
                else:
                    continue # no rows for extvariables.
    return N

def rates_strings(m, fully_qualified=True):
    """Generate a tuple of tuples of
       (name, rate) where
       'name' is the name of a reaction
       'rhs' is the string of the rate of the reaction.
    """
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    return tuple([(v.name, v(fully_qualified=fully_qualified)) for v in m.reactions])

def dXdt_strings(m):
    """Generate a tuple of tuples of
       (name, rhs) where
       'name' is the name of a variable
       'rhs' is the string of the rhs of that variable
       in the SODE for this model.
    """
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    N = genStoichiometryMatrix(m)
    res = []
    for i,name in enumerate(m.varnames):
        dXdtstring = ''
        for j,v in enumerate(m.reactions):
            coef = N[i,j]
            if coef == 0.0: continue
            ratestring = '(%s)'% v(fully_qualified = True)
            if coef == 1.0:
                ratestring = '+'+ratestring
            else:
                ratestring = "%g*%s" % (coef,ratestring)
                if coef > 0.0:
                    ratestring = '%s%s'%('+', ratestring)
            dXdtstring += ratestring
        res.append((name, dXdtstring))
    return tuple(res)


def Jacobian_strings(m, _scale=1.0):
    """Generate a matrix (list of lists) of strings
       to compute the jacobian for this model.
    
       IMPORTANT: sympy module must be installed!"""

    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    try:
        import sympy
    except ImportError:
        print ('ERROR: sympy module must be installed to generate Jacobian strings')
        raise
    dxdtstrings = [d[1] for d in dXdt_strings(m)]
    nvars = len(dxdtstrings)
    vnames = m.varnames
    symbmap, sympysymbs = _gen_canonical_symbmap(m)
    for i in range(nvars):
        dxdtstrings[i] = _replace_exprs2canonical(dxdtstrings[i], symbmap)
    jfuncs = []
    for i in range(nvars):
        jfuncs.append([])
        ids = identifiersInExpr(dxdtstrings[i])
        if len(ids) == 0:
            for j in range(nvars):
                jfuncs[i].append('0.0')
        else:
            for j in range(nvars):
                varsymb = symbmap[vnames[j]]
                res = eval(dxdtstrings[i], None, sympysymbs)
                res = res * _scale
                dres = str(sympy.diff(res, varsymb))
                if dres == '0':
                    dres = '0.0'
                jfuncs[i].append(dres)
    # back to original ids
    for i in range(nvars):
        for j in range(nvars):
            jfuncs[i][j] = _replace_canonical2exprs(jfuncs[i][j], symbmap)
    return jfuncs
        
def dfdp_strings(m, parnames, _scale=1.0):
    """Generate a matrix (list of lists) of strings
       to compute the partial derivatives of rhs of SODE
       with respect to a list of parameters.
       parnames is a list of parameter names.
    
       IMPORTANT: sympy module must be installed!"""

    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    try:
        import sympy
    except ImportError:
        print ('ERROR: sympy module must be installed to generate Jacobian strings')
        raise
    dxdtstrings = [d[1] for d in dXdt_strings(m)]
    nvars = len(dxdtstrings)
    symbmap, sympysymbs = _gen_canonical_symbmap(m)
    for i in range(nvars):
        dxdtstrings[i] = _replace_exprs2canonical(dxdtstrings[i], symbmap)
    npars = len(parnames)
    jfuncs = []
    for i in range(nvars):
        jfuncs.append([])
        ids = identifiersInExpr(dxdtstrings[i])
        if len(ids) == 0:
            for j in range(npars):
                jfuncs[i].append('0.0')
        else:
            for j in range(npars):
                if parnames[j] not in symbmap:
                    dres = '0.0'
                else:
                    varsymb = symbmap[parnames[j]]
                    res = eval(dxdtstrings[i], None, sympysymbs)
                    res = res * _scale
                    dres = str(sympy.diff(res, varsymb))
                    if dres == '0':
                        dres = '0.0'
                jfuncs[i].append(dres)
    # back to original ids
    for i in range(nvars):
        for j in range(npars):
            jfuncs[i][j] = _replace_canonical2exprs(jfuncs[i][j], symbmap)
    return jfuncs


def _gen_canonical_symbmap(m):
    try:
        import sympy
    except ImportError:
        print ('ERROR: sympy module must be installed to generate Jacobian strings')
        raise
    symbmap = {}
    sympysymbs = {}
    symbcounter = 0
    for x in m.varnames:
        symbname = '_symbol_Id%d'% symbcounter
        symbmap[x] = symbname
        sympysymbs[symbname] = sympy.Symbol(symbname)
        symbcounter += 1
    for p in m.parameters:
        symbname = '_symbol_Id%d'% symbcounter 
        symbmap[p.name] = symbname
        sympysymbs[symbname] = sympy.Symbol(symbname)
        symbcounter += 1
    return symbmap, sympysymbs

def _replace_exprs2canonical(s, symbmap):
    for symb in symbmap:
        symbesc = symb.replace('.', '\.')
##         print 'symb =', symb
##         print 'symbesc =', symbesc
##         print 's =', s
##         s = s.replace(symb, symbmap[symb])
        s = re.sub(r"(?<![_.])\b%s\b(?![_.\[])"%symbesc, symbmap[symb], s)
##         print 's =', s
    return s

def _replace_canonical2exprs(s, symbmap):
    for symb in symbmap:
        s = re.sub(r"(?<![.])\b%s\b" % symbmap[symb], symb, s)
##         s = s.replace(symbmap[symb], symb)
    return s

def add_dSdt_to_model(m, pars):
    """Add sensitivity ODEs to model, according to formula:
    
    dS/dt = df/dx * S + df/dp
    
    m is a model object
    pars are a list of parameter names
    """
    
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    try:
        import sympy
    except ImportError:
        print ('ERROR: sympy module must be installed to generate Jacobian strings')
        raise
    #Find pars that are initial values
    init_of = []
    for p in pars:
        if '.' in p:
            tks = p.split('.')
            if tks[1] in m.init:
                init_of.append(tks[1])
            else:
                init_of.append(None)
        else:
            init_of.append(None)

    J = Jacobian_strings(m)
    dfdpstrs = dfdp_strings(m, pars)
    nvars = len(J)
    npars = len(pars)
    
    symbmap, sympysymbs = _gen_canonical_symbmap(m)

    #create symbols for sensitivities
    Snames = []
    Smatrix = []
    for i in range(nvars):
        Smatrix.append([])
        for j in range(npars):
            Sname = "d_%s_d_%s" % (m.varnames[i], pars[j].replace('.','_'))
            sympysymbs[Sname] = sympy.Symbol(str(Sname))
            Smatrix[i].append(Sname)
            Snames.append((m.varnames[i], pars[j], Sname))
    #compute rhs of sensitivities symbolically
    for i in range(nvars):
        vname = m.varnames[i]
        for j in range(npars):
            #compute string for dS/dt
            if init_of[j] is None:
                resstr = dfdpstrs[i][j]
            else:
                resstr = ''
            # matrix multiplication with strings:
            for k in range(nvars):
                resstr = resstr+ "+(%s)*(%s)"%(J[i][k], Smatrix[k][j])
            resstr = _replace_exprs2canonical(resstr, symbmap)
            #make sympy reduce the expression using _symbs dictionary
            res = eval(resstr, None, sympysymbs)
            dres = str(res)
            if dres == '0':
                dres = '0.0'
            dres = _replace_canonical2exprs(dres, symbmap)
            m.set_variable_dXdt(Smatrix[i][j], dres)
            #setattr(m, Smatrix[i][j], variable(dres))
            if init_of[j] is None:
                m.set_init([(Smatrix[i][j], 0.0)])
                #setattr(m.init, Smatrix[i][j], 0.0)
            else:
                if init_of[j] == vname:
                    m.set_init([(Smatrix[i][j], 1.0)])
                    #setattr(m.init, Smatrix[i][j], 1.0)
                else:
                    m.set_init([(Smatrix[i][j], 0.0)])
                    #setattr(m.init, Smatrix[i][j], 0.0)
    return Snames

def _gen_calc_symbmap(m, with_uncertain = False):
    symbmap = {}
    for i, x in enumerate(m.varnames):
        symbname = "variables[%d]"%i
        symbmap[x] = symbname
    for i, x in enumerate(m.input_variables):
        symbname = "input_variables[%d]"%i
        symbmap[x.name] = symbname
    if with_uncertain:
        for i,u in enumerate(m.with_bounds):
            symbname = "m_Parameters[%d]"%i
            symbmap[u.name] = symbname
    for p in m.parameters:
        if p.bounds and with_uncertain:
            continue
        valuestr = "%g"% p
        symbmap[p.name] = valuestr
##     for i, x in enumerate(m._usable_functions):
##         symbname = "l_functions[%d]"%i
##         symbmap[x] = symbname
        
    return symbmap

def rateCalcString(rateString, symbmap):
    return _replace_exprs2canonical(rateString, symbmap)

def compile_rates(m, collection, with_uncertain=False):
    symbmap = _gen_calc_symbmap(m, with_uncertain=with_uncertain)
    ratestrs = [rateCalcString(v(fully_qualified=True), symbmap) for v in collection]
    ratebytecode = [compile(v, '<string>','eval') for v in ratestrs]
    return ratebytecode

def _get_rates_function(m, with_uncertain):
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    
    # get sizes
    n_inputs = len(m.input_variables)
    n_v = len(m.reactions)
    n_transf = len(m.transformations)
    
    # compile input variables
    invarbytecode = compile_rates(m, m.input_variables, with_uncertain=with_uncertain)
    # create array to hold input-variable values's
    input_variables = empty(n_inputs)
    eniv = list(enumerate(invarbytecode))

    # compile reaction rates
    ratebytecode = compile_rates(m, m.reactions, with_uncertain=with_uncertain)
    # create array to hold reaction rate values
    _v_rates = empty(n_v)
    enre = list(enumerate(ratebytecode))
    
    # compile transformations
    tratebytecode = compile_rates(m, m.transformations, with_uncertain=with_uncertain)
    # create array to hold ttransformation values
    _t_rates = empty(n_transf)
    entr = list(enumerate(tratebytecode))
    
    def f(variables, t):
        m_Parameters = m._Model__m_Parameters
        for i, r in eniv:
            input_variables[i] = eval(r, m._usable_functions, locals())
        for i, r in enre:
            _v_rates[i] = eval(r, m._usable_functions, locals())
        for i, r in entr:
            _t_rates[i] = eval(r, m._usable_functions, locals())
        return input_variables, _v_rates, _t_rates
        
    return f
    
def all_rates_func(m, with_uncertain=False, scale=1.0, t0=0.0):
    """Generate function to compute rate vector for this model.
    
       Function has signature f(variables, t)"""
    
    get_rates = _get_rates_function(m, with_uncertain=with_uncertain)
    
    def fout(variables, t):
        t = t*scale + t0
        return get_rates(variables, t)

    return fout

    
def rates_func(m, with_uncertain=False, scale=1.0, t0=0.0):
    """Generate function to compute rate vector for this model.
    
       Function has signature f(variables, t)"""
    
    get_rates = _get_rates_function(m, with_uncertain=with_uncertain)
    
    def fout(variables, t):
        t = t*scale + t0
        _, v, _ = get_rates(variables, t)
        return v

    return fout

def transf_func(m, with_uncertain=False, scale=1.0, t0=0.0):
    """Generate function to compute rate vector for this model.
    
       Function has signature f(variables, t)"""

    get_rates = _get_rates_function(m, with_uncertain=with_uncertain)
    
    def fout(variables, t):
        t = t*scale + t0
        _, _, tvalues = get_rates(variables, t)
        return tvalues
    return fout


def getdXdt(m, with_uncertain=False, scale=1.0, t0=0.0):
    """Generate function to compute rhs of SODE for this model.
    
       Function has signature f(variables, t)
       This is compatible with scipy.integrate.odeint"""

    get_rates = _get_rates_function(m, with_uncertain=with_uncertain)

    # compute stoichiometry matrix, scale and transpose
    N  = genStoichiometryMatrix(m)
    N *= scale
    NT = N.transpose()
    x = empty(len(m.varnames))
    
    def fout(variables, t):
        t = t*scale + t0
        input_variables, v, t_values = get_rates(variables, t)
        dot(v,NT,x)
        return x
    return fout


def getJacobian(m, with_uncertain=False, scale=1.0, t0=0.0):
    """Generate function to compute the jacobian for this model.
    
       Function has signature J(variables, t)
       and returns an nvars x nvars numpy array
       IMPORTANT: sympy module must be installed!"""

    Jstrings = Jacobian_strings(m, _scale = scale)
    nvars = len(Jstrings)
    
    #compile rate laws
    symbmap = _gen_calc_symbmap(m, with_uncertain = with_uncertain)
    ratestrs = [[rateCalcString(col, symbmap) for col in line] for line in Jstrings]
    ratebytecode = [[compile(col, '<string>','eval') for col in line] for line in ratestrs]

    def J(variables, t):
        Jarray = empty((nvars,nvars), float)
        t = t*scale + t0
        for i in range(nvars):
            for j in range(nvars):
                Jarray[i,j] = eval(ratebytecode[i][j], m._usable_functions, locals())
        return Jarray
    def J2(variables, t):
        m_Parameters = m._Model__m_Parameters
        Jarray = empty((nvars,nvars), float)
        t = t*scale + t0
        for i in range(nvars):
            for j in range(nvars):
                Jarray[i,j] = eval(ratebytecode[i][j], m._usable_functions, locals())
        return Jarray
    if with_uncertain:
        return J2
    else:
        return J


def _gen_outputs_decl(m, ignore_replist=False):
    decl = m.metadata.get('!!', None)
    
    if decl is not None and not ignore_replist:
        names = decl.strip().split()
        return genTransformationFunction(m, names)
    else:
        return None

def genTransformationFunction(m, f):
    special_transf = ['~']
    special_rates = ['>', '>>', '->']
    all_special = special_transf + special_rates
    
    if f in all_special:
        f = [f.strip()]
    
    if _is_sequence(f):
        names = []
        for a in f:
            if not _is_string(a):
                raise TypeError(str(a) + ' must be a string')
            if a in special_transf:
                names.extend([x.name for x in m.transformations])
            elif a in special_rates:
                names.extend([x.name for x in m.reactions])
            else:
                names.append(a)
    else:
        raise TypeError('outputs must be a sequence of names.')
    
    nargs = len(names)
    
    get_rates = _get_rates_function(m, with_uncertain=False)

    data = []

    for name in names:
        i, kind = _find_indexof_component(m, name)
        if kind == 'parameters':
            data.append(('p', m.getp(name)))
        elif kind == 'variables':
            data.append(('v', i))
        elif kind == 'transf':
            data.append(('t', i))
        elif kind == 'reactions':
            data.append(('r', i))
        elif kind == 'invar':
            data.append(('i', i))
        elif kind == 'uncertain':
            data.append(('u', i))
        else:
            raise AttributeError(name + ' is not a component in this model')
    
    args = [0.0] * nargs
    for (i, (kind, d)) in enumerate(data):
        if kind == 'p':
            args[i] = d

    def fout(variables, t):
        input_variables, v, t_values = get_rates(variables, t)
        
        for (i, (kind, d)) in enumerate(data):
            if kind == 'p':
                continue
            elif kind == 'v':
                args[i] = variables[d]
            elif kind == 'r':
                args[i] = v[d]
            elif kind == 'i':
                args[i] = input_variables[d]
            elif kind == 't':
                args[i] = t_values[d]
            elif kind == 'u':
                args[i] = m._Model__m_Parameters[d]
            else:
                continue
        return args
            
    result = fout
    result.names = names
    return result


def solve(model, 
          tf=None, 
          npoints=500, 
          t0=0.0, 
          initial='init', 
          times=None, 
          outputs=None, 
          title=None,
          ignore_replist=False):
    
    solver=integrate._odepack.odeint
    names = [x for x in model.varnames]

    #get initial values
    if initial == 'init':
        y0 = init2array(model)
    else:
        y0 = copy(initial)
    if tf is None:
        try:
            tf = float(model.metadata.get('tf', None))
        except:
            pass
        if tf is None:
            tf = 1.0
    if times is None:
        times = linspace(t0, tf, npoints)            

    # scale times to maximum time in data
    t0 = times[0]
    scale = float(times[-1] - t0)
    #scale = 1.0
    
    f = getdXdt(model, scale=scale, t0=t0)
    t = copy((times-t0)/scale)  # this scales time points
    
    output = solver(f, y0, t, (), None, 0, -1, -1, 0, None, 
                    None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
    if output[-1] < 0: return None
    Y = output[0]
    
    if title is None:
        title = model.metadata.get('title', '')        
    Y = copy(Y.T)

    sol = SolutionTimeCourse (times, Y, names, title, dense = True)
    
    # get outputs
    f = _gen_outputs_decl(model, ignore_replist)
    # overide if outputs argument is not None
    if outputs is not None:
        f = genTransformationFunction(model, outputs)
    if f is not None:
        sol.apply_transf(f, f.names)
    
    return sol

class ModelSolver(object):
    def __init__(self,
          model, 
          tf=1.0, 
          npoints=500, 
          t0=0.0, 
          initial='init', 
          times=None, 
          outputs=None, 
          title=None,
          ignore_replist=False,
          changing_pars=None):
        
        self.model = model.copy()
        # reset all bounds
        bnames = [p.name for p in self.model.with_bounds]
        for name in bnames:
            self.model.reset_bounds(name)
        
        self.names = [x for x in self.model.varnames]
        self.title = title
        if self.title is None:
            self.title = self.model.metadata.get('title', '')
        
        #get initial values
        if initial == 'init':
            self.y0 = copy(init2array(self.model))
        else:
            self.y0 = copy(initial)
        
        self.times = times
        if self.times is None:
            self.times = linspace (t0, tf, npoints)
        
        # scale times to maximum time in data
        t0 = self.times[0]
        scale = float(self.times[-1] - t0)
        self.t  = (self.times-t0)/scale  # this scales time points

        # store names of changing parameters
        if changing_pars is None:
            changing_pars = []
        if _is_string(changing_pars):
            changing_pars = changing_pars.strip().split()
        self.changing_pars = changing_pars
        
        # find initial values in changing parameters
        mapinit2pars = []
        for i, parname in enumerate(self.changing_pars):
            if parname.startswith('init'):
                varname = parname.split('.')[-1]
                ix = self.model.varnames.index(varname)
                mapinit2pars.append((ix,i))
            self.model.set_bounds(parname, (0,1)) # bogus bounds
                
        self.pars_initindexes = array([j for (i,j) in mapinit2pars], dtype=int)
        self.vars_initindexes = array([i for (i,j) in mapinit2pars], dtype=int)

        self.f = getdXdt(self.model, with_uncertain=True, scale=scale, t0=t0)
        
        self.tranf_f = None
        self.tranf_names = None
        
        # get outputs
        f = _gen_outputs_decl(self.model, ignore_replist)
        # overide if outputs argument is not None
        if outputs is not None:
            f = genTransformationFunction(self.model, outputs)
        if f is not None:
            self.tranf_f = f
            self.tranf_names = f.names

    def solve(self, title = None, par_values = None):
        
        # set initial values
        y0 = copy(self.y0)
        
        # set varying parameters (may be initial values)
        if par_values is not None:
            par_values = array(par_values)
            self.model.set_uncertain(par_values)
            # fill uncertain initial values
            y0[self.vars_initindexes] = par_values[self.pars_initindexes]
                
        
        output = integrate._odepack.odeint(self.f, y0, self.t, (), 
                                           None, 0, -1, -1, 0, None,
                                           None, None, 0.0, 0.0, 0.0, 
                                           0, 0, 0, 12, 5)
        if output[-1] < 0: return None
        Y = output[0]
        if title is None:
            title = self.title
        sol = SolutionTimeCourse(self.times, Y.T, self.names, title, dense=True)
        
        # a filter string or transformation function
        if self.tranf_f is not None:
            sol.apply_transf(self.tranf_f, self.tranf_names)
        return sol

def scan(model, plan,
        tf = 1.0, 
        npoints = 500, 
        t0 = 0.0, 
        initial = 'init', 
        times=None, 
        outputs=None, 
        titles=None,
        changing_pars = None):
        
        """Wrapper around ModelSolver."""
                        
        plan = dict(plan)
        names = list(plan.keys())
        # zip, terminating on the shortestsequence
        design = list(izip(*list(plan.values())))
        nruns = len(design)
        
        if titles is None:
            titles = []
            for v in design:
                pairs = list(zip(names, v))
                pairs = ['%s = %g'%(n,v) for (n,v) in pairs]
                titles.append(', '.join(pairs))

##         print 'plan'
##         print plan
##         print '---'
##         print 'parameter names'
##         print names
##         print '---'
##         print 'design'
##         print design
##         print '---'
##         print 'titles'
##         print titles
##         print '---'
        
        ms = ModelSolver(model, tf=tf, npoints=npoints, t0=t0, 
                        initial=initial, times=times, outputs=outputs, 
                        changing_pars=names)
    
        s = Solutions()
        for title,p in zip(titles,design):
            s += ms.solve(title = title, par_values = p)

        return s
    

def test():
    #import time
    from modelparser import read_model     
    m = read_model("""
    title a simple 2 enzyme system
    v1 = A -> B, rate = V*A/(Km1 + A), V = 1, Km = 1
    v2 = B ->  , rate = V*B/(Km2 + B)
    V  = sqrt(4.0)
    Km1 = 1
    Km2 = 0.2
    find Km2 in [0, 1.2]
    init = state(B = 0.4, A = 1)
    -> vin = 2 * A * v1.Km
    ~ t1 = A+B + vin
    ~ t2 = v1.V * A * step(t, 1.0)
    # ~ t3 = v1.V * A * max(t, 1.0)
    """)
    print (m)

    print ('********** Testing stoichiometry matrix ********************')
    print ('Stoichiometry matrix:')
    N = genStoichiometryMatrix(m)
    print ('  ', '  '.join([v.name for v in m.reactions]))
    for i,x in enumerate(m.varnames):
        print (x, N[i, :])
    print()
    print ('********** Testing state2array()****************************')
    print ('state2array(m):')
    v = init2array(m)
    print (v, 'of type', type(v))
    print()
    print ('********** Testing rate and dXdt strings *******************')
    print ('rates_strings(fully_qualified = False): ---')
    for v in rates_strings(m, fully_qualified = False):
        print (v)
    print ('\nrates_strings(): -------------------------')
    for v in rates_strings(m):
        print (v)
    print ('\ndXdt_strings(): --------------------------')
    for xname,dxdt in dXdt_strings(m):
        print ('(d%s/dt) ='%(xname),dxdt)
    print()
    print ('Jacobian_strings(): -------------------------')
    vnames = m.varnames
    for i,vec in enumerate(Jacobian_strings(m)):
        for j, dxdx in enumerate(vec):
            print ('(d d%s/dt / d %s) ='%(vnames[i],vnames[j]), dxdx)
    print()
    print ('dfdp_strings(m, parnames): ------------------')
    parnames = "Km2 v1.V".split()
    print ('parnames = ["Km2", "v1.V"]\n')
    vnames = m.varnames
    for i,vec in enumerate(dfdp_strings(m, parnames)):
        for j, dxdx in enumerate(vec):
            print ('(d d%s/dt / d %s) ='%(vnames[i],parnames[j]), dxdx)
    print()
    print ('dfdp_strings(m, parnames): (with unknown pars)')
    parnames = "Km3 v1.V".split()
    print ('parnames = ["Km3", "v1.V"]\n')
    vnames = m.varnames
    for i,vec in enumerate(dfdp_strings(m, parnames)):
        for j, dxdx in enumerate(vec):
            print ('(d d%s/dt / d %s) ='%(vnames[i],parnames[j]), dxdx)
    print ('\n********** Testing _gen_calc_symbmap(m) *******************')
    print ('_gen_calc_symbmap(m, with_uncertain = False):')
    print (_gen_calc_symbmap(m))
    print ('\n_gen_calc_symbmap(m, with_uncertain = True):')
    print (_gen_calc_symbmap(m, with_uncertain = True))
    
    print ('\n********** Testing rateCalcString **************************')
    symbmap = _gen_calc_symbmap(m, with_uncertain = False)
    symbmap2 = _gen_calc_symbmap(m, with_uncertain = True)
    for v in (m.reactions.v1, 
              m.reactions.v2, 
              m.transformations.t1, 
              m.transformations.t2,
              m.input_variables.vin):
        vstr = v(fully_qualified = True)
        print ('calcstring for %s = %s\n\t'% (v.name, vstr), rateCalcString(vstr, symbmap))
    print ('calcstring for v2 with uncertain parameters:\n\t', rateCalcString(m.reactions.v2(fully_qualified = True), symbmap2))

    print ('\n********** Testing rate and dXdt generating functions ******')
    print ('Operating point --------------------------------')
    varvalues = [1.0, 0.4]
    pars      = [0.4]
    t         = 0.0
    
    dxdtstrs = [b for (a,b) in dXdt_strings(m)]

    print ("t =", t)
    print ('variables:')
    print (dict((n, value) for n,value in zip(m.varnames, varvalues)))
    print ('parameters:')
    print (dict((p.name, p) for p in m.parameters))
 
    print ('\n---- rates using rates_func(m) -------------------------')
    vratesfunc = rates_func(m)
    vrates = vratesfunc(varvalues,t)
    frmtstr = "%s = %-25s = %s"
    for v,r in zip(m.reactions, vrates):
        print (frmtstr % (v.name, v(fully_qualified = True), r))

    print ('---- transformations using transf_func(m) ----------------')
    tratesfunc = transf_func(m)
    trates = tratesfunc(varvalues,t)
    for v,r in zip(m.transformations, trates):
        print (frmtstr % (v.name, v(fully_qualified = True), r))
    print ('---- same, at t = 2.0 --')
    trates = tratesfunc(varvalues,2.0)
    for v,r in zip(m.transformations, trates):
        print (frmtstr % (v.name, v(fully_qualified = True), r))

    print ('--------- NEW MODEL, with input variables -------------')
    
    m2 = read_model("""
    title a simple 2 enzyme system
    v1 = A -> B, rate = vin*A/(Km1 + A), V = 1, Km = 1
    v2 = B ->  , rate = V*B/(Km2 + B)
    V  = sqrt(4.0)
    Km1 = 1
    Km2 = 0.2
    find Km2 in [0, 1.2]
    init = state(B = 0.4, A = 1)
    
    -> vin = v1.V * step(t, 1.0) + 1
    """)
    
    print (m2)

    print ('---- dXdt using getdXdt(m) -------------------')
    f = getdXdt(m2)
    dXdt = f(varvalues,t)
    for x,r in zip(m2.varnames,  dXdt):
        print ("d%s/dt = %s" % (x,r))

    print ('---- dXdt using getdXdt(m) setting uncertain parameters ---')
    print ('f = getdXdt(m, with_uncertain = True)')
    f = getdXdt(m2, with_uncertain = True)
    print ('setting uncertain as', dict((v.name, value) for v,value in zip(m2.with_bounds, pars)))
    print ('m.set_uncertain(pars)')
    m2.set_uncertain(pars)
    dXdt = f(varvalues,t)
    for x,r in zip(m2.varnames, dXdt):
        print ("d%s/dt = %s" % (x, r))

    print ('---- dXdt using getdXdt(m) with a state argument (m.init) --')
    print ('m.init:', m2.get_init())
    print ('f = getdXdt(m)')
    f = getdXdt(m2)
    print ('dXdt = f(init2array(m),t)')
    dXdt = f(init2array(m2),t)
    for x,r in zip(m2.varnames, dXdt):
        print ("d%s/dt = %s" % (x, r))
    
    print ('---- same, changing state argument ---------------------------')
    m2.init.A = 2.0
    print ('after m2.init.A = 2.0')
    print ('m2.init:', m2.get_init())
    print ('dXdt = f(init2array(m2),t)')
    dXdt = f(init2array(m2),t)
    for x,r in zip(m2.varnames, dXdt):
        print ("d%s/dt = %s" % (x, r))
    print ('\n********** Testing add_dSdt_to_model functions ***************')
    print ('\n--- Using back model m ---\n')
    m2 = m.copy()
    print (m2)
    print ('----------------------------------------------------')
    pars = "Km2 v1.V".split()
    print ('pars =', pars)
    print ("\n!!! applying function add_dSdt_to_model(m, pars) !!!")
    Snames = add_dSdt_to_model(m2, pars)
    print ('Snames = \n', Snames)
    print (m2)
    
    # print 'BEGIN EXAMPLES'
    # t0 = time.time()
    
    print ('---------------- EXAMPLE 1 ------------------')
    mtext = """
    title a simple 2 enzyme system
    v1 : A -> B, rate = Vin*A/(Km + A), V = 0.1, Km = 1
    v2 : B -> C, rate = V*B/(Km + B), V = sqrt(4.0), Km = 20

    init : A = 1
    ~ sum = A + B + C
    ~ sumAB = A + B
    -> Vin = 0.1 * step(t, 10)
    !! A B C ~
    """
    print (mtext)

    m1 = read_model(mtext)

    solution1 = solve(m1, tf=50, title='two enzymes, use !! A C ~')
    solution1a = solve(m1, tf=50, outputs='A B C sum'.split(),
                       title='explicit outputs=[A B C sum]')
    solution1v = solve(m1, tf=100, outputs='>>',
                       title='outputs=">>"')

    print ('--- Last time point ----')
    print ('At t =', solution1.t[-1])
    #print solution1.last
    for x in solution1.last:
        print ("%-8s= %f" % (x, solution1.last[x]))
    
    # print 'END of EXAMPLES 1'
    # t1 = time.time()
    # print 'took', t1 - t0
    
    print ('---------------- EXAMPLE 3 ------------------')
    m3 = read_model(models.ca.text)

    print (models.ca.text)
    ms = ModelSolver(m3, tf = 8.0, npoints = 2000)
    solution3 = ms.solve()
##     solution3 = solve(m3, tf = 8.0, npoints = 2000)

    # print 'END of EXAMPLES 3'
    # t3 = time.time()
    # print 'took', t3 - t1


    print ('---------------- EXAMPLE 4 ------------------')
    m4 = read_model(models.rossler.text)

    print (m4)

    solution4 = solve(m4, tf = 100.0, npoints = 2000, 
                      outputs="x1 x2 x3".split())
    solution4b = solve(m4, tf = 100.0, npoints = 2000, outputs="~",
                      title='Rossler, outputs="~"')
    
    def transformation(vars,t):
        if t > 40.0:
            return (vars[0]-5.0, vars[1], vars[2])
        else:
            return (-5.0, vars[1], vars[2])

    solution4.apply_transf(transformation,
                           new_title='Rossler, after a transformation')

    # print 'END of EXAMPLES 4'
    # t4 = time.time()
    # print 'took', t4 - t3
    
    #savingfile = open('examples/analysis.png', 'w+b')
    savingfile = 'examples/analysis.png'
    sols = Solutions([solution1, solution1a, solution1v,
                      solution3, 
                      solution4b, solution4])
    sols.plot() #save2file=savingfile)
    
    # print 'END of plotting first 4 examples'
    # tplot = time.time()
    # print 'took', tplot - t4


    print ('---------------- scanning example ------------------')
    m3 = read_model(models.ca.text)
    scans = 0.0, 0.1, 0.3, 0.5, 0.8, 1.0
    
    sols2 = scan(m3, {'B': scans}, tf=10.0)
    # print 'END of SCANNING EXAMPLE'
    # tscancomp = time.time()
    # print 'took', tscancomp - tplot

    sols2.plot(legend=True, ynormalize=True,  group=['Ca'],
               fig_size=(16,9))

    # print 'END of PLOTTING SCANNING EXAMPLE'
    # tscan = time.time()
    # print 'took', tscan - tscancomp



    print ('---------------- stairway example ------------------')
    mtext = """
    title a simple 2 enzyme system
    v1 : A -> B, rate = Vin*A/(Km + A), V = 0.1, Km = 1
    v2 : B -> C, rate = V*B/(Km + B), V = 10, Km = 20
    v3 : C ->, rate = kout * C, kout = 1
    A = 1

    init : B = 0, C = 0

    -> Vin = stairway(t, [50, 100, 150, 200, 250], [1, 2, 3, 4, 5])
    !! Vin B C
    """
    print (mtext)

    mstair = read_model(mtext)

    solstairs = solve(mstair, tf=300, title='stairway')
    # print 'END of STAIRWAY EXAMPLE'
    # tstairway = time.time()
    # print 'took', tstairway - tscan

    solstairs.plot(fig_size=(16,9), show=True)

    # print 'END of STAIRWAY PLOTTING'
    # tstairwayplot = time.time()
    # print 'took', tstairwayplot - tstairway


if __name__ == "__main__":
    test()
 