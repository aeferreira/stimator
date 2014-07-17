#!/usr/bin/env python
# -*- coding: utf-8- *-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator dynamical systems related functions
# Copyright António Ferreira 2006-2014
#----------------------------------------------------------------------------
from model import *
import math
from kinetics import *
from numpy import *
import pprint
from modelparser import read_model
import os
import os.path
from scipy import integrate
from timecourse import SolutionTimeCourse, Solutions

from examples import models


def state2array(m, state):
    """Transforms a state object into a numpy.array object.
       
       This is necessary for most numerical functions of numpy+scipy.
       Can accept the name of a state (must exist in Model) or state object.
       Values are returned in the order of model variables.
    """
    if isinstance(state, str) or isinstance(state, unicode):
        if not hasattr(m, state):
            raise AttributeError( state + ' is not defined for this model')
        state = getattr(m, state)
    return array([state._ownparameters.get(var,0.0) for var in m().varnames])

def genStoichiometryMatrix(m):
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)

    vnames = m().varnames
    N = zeros((len(vnames),len(m().reactions)), dtype=float)
    for j,v in enumerate(m().reactions):
        for rORp, signedunit in [(v._reagents,-1.0),(v._products,1.0)]:
            for c in rORp:
                coef, var = (c[1]*signedunit, c[0])
                if var in vnames:
                    ivariable = vnames.index(var) # error handling here
                    N[ivariable, j] = coef
                else:
                    continue # there are no rows for extvariables in stoich. matrix
    return N

def rates_strings(m, fully_qualified = True):
    """Generate a tuple of tuples of
       (name, rate) where
       'name' is the name of a reaction
       'rhs' is the string of the rate of the reaction.
    """
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    return tuple([(get_name(v), v(fully_qualified = fully_qualified)) for v in m().reactions])

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
    for i,name in enumerate(m().varnames):
        dXdtstring = ''
        for j,v in enumerate(m().reactions):
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

def _gen_canonical_symbmap(m):
    try:
        import sympy
    except:
        print 'ERROR: sympy module must be installed to generate Jacobian strings'
        raise
    symbmap = {}
    sympysymbs = {}
    symbcounter = 0
    for x in m().varnames:
        symbname = '_symbol_Id%d'% symbcounter
        symbmap[x] = symbname
        sympysymbs[symbname] = sympy.Symbol(symbname)
        symbcounter += 1
    for p in m().parameters:
        symbname = '_symbol_Id%d'% symbcounter 
        symbmap[get_name(p)] = symbname
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
        s = re.sub(r"(?<![.])\b%s\b"%symbmap[symb], symb, s)
##         s = s.replace(symbmap[symb], symb)
    return s

def Jacobian_strings(m, _scale = 1.0):
    """Generate a matrix (list of lists) of strings
       to compute the jacobian for this model.
    
       IMPORTANT: sympy module must be installed!"""

    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    try:
        import sympy
    except:
        print 'ERROR: sympy module must be installed to generate Jacobian strings'
        raise
    dxdtstrings = [d[1] for d in dXdt_strings(m)]
    nvars = len(dxdtstrings)
    vnames = m().varnames
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
        
def dfdp_strings(m, parnames, _scale = 1.0):
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
    except:
        print 'ERROR: sympy module must be installed to generate Jacobian strings'
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
    except:
        print 'ERROR: sympy module must be installed to generate Jacobian strings'
        raise
    #Find pars that are initial values
    init_of = []
    for p in pars:
        if '.' in p:
            tks = p.split('.')
            if isinstance(getattr(m,tks[0]), StateArray):
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
            Sname = "d_%s_d_%s" % (m().varnames[i], pars[j].replace('.','_'))
            sympysymbs[Sname] = sympy.Symbol(str(Sname))
            Smatrix[i].append(Sname)
            Snames.append((m().varnames[i], pars[j], Sname))
    #compute rhs of sensitivities symbolically
    for i in range(nvars):
        vname = m().varnames[i]
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
            setattr(m, Smatrix[i][j], variable(dres))
            if init_of[j] is None:
                setattr(m.init, Smatrix[i][j], 0.0)
            else:
                if init_of[j] == vname:
                    setattr(m.init, Smatrix[i][j], 1.0)
                else:
                    setattr(m.init, Smatrix[i][j], 0.0)
    return Snames

def _gen_calc_symbmap(m, with_uncertain = False):
    symbmap = {}
    for i, x in enumerate(m().varnames):
        symbname = "variables[%d]"%i
        symbmap[x] = symbname
    if with_uncertain:
        for i,u in enumerate(m().uncertain):
            symbname = "m_Parameters[%d]"%i
            symbmap[get_name(u)] = symbname
    for p in m().parameters:
        if p.bounds and with_uncertain:
            continue
        symbname = "%g"% p
        symbmap[get_name(p)] = symbname
    return symbmap

def rateCalcString(rateString, symbmap):
    return _replace_exprs2canonical(rateString, symbmap)

def compile_rates(m, collection, with_uncertain = False):
    symbmap = _gen_calc_symbmap(m, with_uncertain = with_uncertain)
    ratestrs = [rateCalcString(v(fully_qualified = True), symbmap) for v in collection]
    ratebytecode = [compile(v, '<string>','eval') for v in ratestrs]
    return ratebytecode
    
def rates_func(m, with_uncertain = False, transf = False, scale = 1.0, t0=0.0):
    """Generate function to compute rate vector for this model.
    
       Function has signature f(variables, t)"""

    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    if transf :
        collection = m().transformations
    else:
        collection = m().reactions

    #compile rate laws
    ratebytecode = compile_rates(m, collection, with_uncertain = with_uncertain)
    # create array to hold v's
    v = empty(len(collection))
    en = list(enumerate(ratebytecode))
        
    def f(variables, t):
        t = t*scale + t0
        for i,r in en:
            v[i] = eval(r)
        return v
    def f2(variables, t):
        m_Parameters = m._Model__m_Parameters
        t = t*scale + t0
        for i,r in en:
            v[i] = eval(r)
        return v

    if with_uncertain:
        return f2
    else:
        return f

def genTransformationFunction(m, f):
    if not callable(f):
        if isinstance(f,list) or isinstance(f, tuple):
            for a in f:
                if not (isinstance(a,str) or isinstance(a,unicode)):
                    raise TypeError(str(a) + ' must be a list of strings')
            argnames = f[:]
            names = argnames
        else:
            if not (isinstance(f,str) or isinstance(f,unicode)):
                raise TypeError('argument must be a string, list or callable.')
            argnames = f.split()
            names = argnames
        nargs = len(argnames)
    else:
        cc = f.func_code
        nargs = cc.co_argcount
        argnames = cc.co_varnames[:nargs]
        names = list(argnames[:])
        if hasattr(f, 'names'):
            names[:len(f.names)] = f.names
    data = []
    symbmap = _gen_calc_symbmap(m, with_uncertain = False)
    for a in argnames:
        i, collection = m._Model__findComponent(a)
        if collection == 'parameters':
            data.append(('p', getattr(m, a)))
        elif collection == 'variables':
            data.append(('v', i))
        elif collection == 'transf':
            vstr = rateCalcString(m().transformations[i](fully_qualified = True), symbmap)
            data.append(('t', compile(vstr, '<string>','eval')))
        elif collection == 'reactions':
            vstr = rateCalcString(m().reactions[i](fully_qualified = True), symbmap)
            data.append(('r', compile(vstr, '<string>','eval')))
        else:
            raise AttributeError(a + ' is not a component in this model')
    args = [0.0]*nargs
    en = [(i,c,d) for (i, (c,d)) in enumerate(data)]
    for i,c,d in en:
        if c == 'p':
            args[i] =d

    def retf(variables, t):
        for i,c,d in en:
            if c == 'p':
                continue
            elif c == 'v':
                args[i] = variables[d]
            else:
                args[i] = eval(d)
        return f(*args)
    def retargs(variables, t):
        for i,c,d in en:
            if c == 'p':
                continue
            elif c == 'v':
                args[i] = variables[d]
            else:
                args[i] = eval(d)
        return args
            
    if callable(f):
        result = retf
        result.names = names
        return result
    else:
        result = retargs
        result.names = names
        return result

def getdXdt(m, with_uncertain = False, scale = 1.0, t0=0.0):
    """Generate function to compute rhs of SODE for this model.
    
       Function has signature f(variables, t)
       This is compatible with scipy.integrate.odeint"""

    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    #compile rate laws
    ratebytecode = compile_rates(m, m().reactions, with_uncertain = with_uncertain)

    # compute stoichiometry matrix, scale and transpose
    N  = genStoichiometryMatrix(m)
    N *= scale
    NT = N.transpose()
    # create array to hold v's
    v = empty(len(m().reactions))
    x = empty(len(m().varnames))
    en = list(enumerate(ratebytecode))
    
    def f2(variables, t):
        m_Parameters = m._Model__m_Parameters
        t = t*scale + t0
        for i,r in en:
            v[i] = eval(r)
        dot(v,NT,x)
        return x
    return f2

def getdXdt_exposing_rbc(m, expose_enum, with_uncertain = False, scale = 1.0, t0=0.0, changing_pars = None):
    """Generate function to compute rhs of SODE for this model.
    
       Function has signature f(variables, t)
       This is compatible with scipy.integrate.odeint"""

    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    #compile rate laws
    ratebytecode = compile_rates(m, m().reactions, with_uncertain = with_uncertain)
    # compute stoichiometry matrix, scale and transpose
    N  = genStoichiometryMatrix(m)
    N *= scale
    NT = N.transpose()
    # create array to hold v's
    v = empty(len(m().reactions))
    x = empty(len(m().varnames))
    for i in range(len(m().reactions)):
        expose_enum[i] = (i,ratebytecode[i])
    
    def f2(variables, t):
        m_Parameters = m._Model__m_Parameters
        t = t*scale + t0
        for i,r in expose_enum:
            v[i] = eval(r)
        dot(v,NT,x)
        return x
    return f2


def dXdt_with(m, uncertainparameters, scale = 1.0, t0=0.0):
    """Generate function to compute rhs of SODE for this model.
    
       Function has signature f(variables, t)
       This is compatible with scipy.integrate.odeint"""

    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    #compile rate laws
    ratebytecode = compile_rates(m, m().reactions, with_uncertain = True)
    # compute stoichiometry matrix, scale and transpose
    N  = genStoichiometryMatrix(m)
    N *= scale
    NT = N.transpose()

    # create array to hold v's
    v = empty(len(m().reactions))
    en = list(enumerate(ratebytecode))
    def f(variables, t):
        m_Parameters = uncertainparameters
        t = t*scale + t0
        for i,r in en:
            v[i] = eval(r)
        return dot(v,NT)
    return f

def getJacobian(m, with_uncertain = False, scale = 1.0, t0=0.0):
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
                Jarray[i,j] = eval(ratebytecode[i][j])
        return Jarray
    def J2(variables, t):
        m_Parameters = m._Model__m_Parameters
        Jarray = empty((nvars,nvars), float)
        t = t*scale + t0
        for i in range(nvars):
            for j in range(nvars):
                Jarray[i,j] = eval(ratebytecode[i][j])
        return Jarray
    if with_uncertain:
        return J2
    else:
        return J

def solve(model, 
          tf = None, 
          npoints = 500, 
          t0 = 0.0, 
          initial = 'init', 
          times = None, 
          outputs=False, 
          title = None,
          ignore_replist = False):
    
    salg=integrate._odepack.odeint
    names = [x for x in model().varnames]

    #get initial values, possibly from a state in the model
    if isinstance(initial, str) or isinstance(initial, StateArray):
        y0 = copy(state2array(model,initial))
    else:
        y0 = copy(initial)
    if tf is None:
        try:
            tf = float(model['tf'])
        except:
            pass
        if tf is None:
            tf = 1.0
    if times is None:
        times = linspace (t0, tf, npoints)
##     elif 
##         if isinstance(times, list):
            
    
    # scale times to maximum time in data
    t0 = times[0]
    scale = float(times[-1] - t0)
    #scale = 1.0
    
    f = getdXdt(model, scale=scale, t0=t0)
    t = copy((times-t0)/scale)  # this scales time points
    output = salg(f, y0, t, (), None, 0, -1, -1, 0, None, 
                    None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
    if output[-1] < 0: return None
    Y = output[0]
    if title is None:
        title = model['title']        
    Y = copy(Y.T)

    sol = SolutionTimeCourse (times, Y, names, title)
    
    #get outputs
    if outputs is False: # variables are output
        pass
    elif outputs is True: #transformations are output
        #compute model transformations
        f     = model.transf_func()
        names = [get_name(x) for x in model.transf]
        sol.apply_transf(f,names)
    else:
        #a filter string or transformation function
        f = genTransformationFunction(model, outputs)
        sol.apply_transf(f, f.names)
    if model['!!'] is not None and ignore_replist == False:
        names = model['!!'].split()
        sol = sol.copy(names=names)
    
    return sol

class ModelSolver(object):
    def __init__(self,
          model, 
          tf = 1.0, 
          npoints = 500, 
          t0 = 0.0, 
          initial = 'init', 
          times = None, 
          outputs=False, 
          title = None,
          changing_pars = None):
        self.model = model.clone()
        self.salg=integrate._odepack.odeint
        self.names = [x for x in self.model().varnames]
        self.changing_pars = changing_pars
        if self.changing_pars is None:
            self.changing_pars = []
        if isinstance(self.changing_pars, str):
            self.changing_pars = self.changing_pars.split()
        self.par_enumeration = []

        # find initial values in changing parameters
        mapinit2pars = []
        for ipar, parname in enumerate(self.changing_pars):
            if parname.startswith('init'):
                varname = parname.split('.')[-1]
                ix = self.model().varnames.index(varname)
                mapinit2pars.append((ix,ipar))
            else:
                self.par_enumeration.append((ipar,parname))
                
        self.pars_initindexes = array([j for (i,j) in mapinit2pars], dtype=int)
        self.vars_initindexes = array([i for (i,j) in mapinit2pars], dtype=int)

        #get initial values, possibly from a state in the model
        if isinstance(initial, str) or isinstance(initial, StateArray):
            self.y0 = copy(state2array(self.model,initial))
        else:
            self.y0 = copy(initial)
        self.times = times
        if self.times is None:
            self.times = linspace (t0, tf, npoints)
        
        # scale times to maximum time in data
        t0 = self.times[0]
        self.t0 = t0
        scale = float(self.times[-1] - t0)
        self.scale = scale
        #scale = 1.0
        self.expose_enum = [0.0 for i in range(len(self.model().reactions))]
        self.f = getdXdt_exposing_rbc(self.model, self.expose_enum, scale=scale, t0=t0)
        self.t  = (self.times-t0)/scale  # this scales time points
        self.title = title
        if self.title is None:
            self.title = self.model['title']
        self.tranf_f = None
        self.tranf_names = None
        #get outputs
        if outputs is False: # variables are output
            pass
        elif outputs is True: #transformations are output
            #compute model transformations
            self.tranf_f     = self.model.transf_func()
            self.tranf_names = [get_name(x) for x in self.model.transf]
        elif isinstance(outputs, str) or callable(outputs): 
            #a filter string or transformation function
            self.tranf_f = genTransformationFunction(self.model, outputs)
            self.tranf_names = self.tranf_f.names
        else:
            raise TypeError("'outputs' parameter is of the wrong type")
    
    def solve(self, title = None, par_values = None):
        y0 = copy(self.y0)
        if par_values is not None:
            par_values = array(par_values)
            for (ip, pname) in self.par_enumeration:
                setattr(self.model, pname, par_values[ip])
            y0[self.vars_initindexes] = par_values[self.pars_initindexes]
        vs = self.model().reactions
        ratebytecode = compile_rates(self.model, vs, with_uncertain = False)
        for i, rbc in enumerate(ratebytecode):
            self.expose_enum[i] = (i,rbc)
        output = integrate._odepack.odeint(self.f, y0, self.t, (), None, 0, -1, -1, 0, None, 
                        None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
        if output[-1] < 0: return None
        Y = output[0]
        if title is None:
            title = self.title
        sol = SolutionTimeCourse (self.times, Y.T, self.names, title)
        
        if self.tranf_f is not None:
            sol.apply_transf(self.tranf_f, self.tranf_names)
        return sol

def test():
     
    m = read_model("""
    title a simple 2 enzyme system
    v1 = A -> B, rate = V*A/(Km1 + A), V = 1, Km = 1
    v2 = B ->  , rate = V*B/(Km2 + B)
    V  = sqrt(4.0)
    Km1 = 1
    Km2 = 0.2
    find Km2 in [0, 1.2]
    init = state(B = 0.4, A = 1)
    ~ t1 = A+B
    ~ t2 = v1.V * A * step(t, 1.0)
    """)
    print m

    print '********** Testing stoichiometry matrix ********************'
    print 'Stoichiometry matrix:'
    N = genStoichiometryMatrix(m)
    print '  ', '  '.join([get_name(v) for v in m().reactions])
    for i,x in enumerate(m().varnames):
        print x, N[i, :]
    print
    print '********** Testing state2array()****************************'
    print 'state2array(m,"init"):'
    v = state2array(m,"init")
    print v, 'of type', type(v)
    print
    print '********** Testing rate and dXdt strings *******************'
    print 'rates_strings(fully_qualified = False): ---'
    for v in rates_strings(m, fully_qualified = False):
        print v
    print '\nrates_strings(): -------------------------'
    for v in rates_strings(m):
        print v
    print '\ndXdt_strings(): --------------------------'
    for xname,dxdt in dXdt_strings(m):
        print '(d%s/dt) ='%(xname),dxdt
    print
    print 'Jacobian_strings(): -------------------------'
    vnames = m().varnames
    for i,vec in enumerate(Jacobian_strings(m)):
        for j, dxdx in enumerate(vec):
            print '(d d%s/dt / d %s) ='%(vnames[i],vnames[j]), dxdx
    print
    print 'dfdp_strings(m, parnames): ------------------'
    parnames = "Km2 v1.V".split()
    print 'parnames = ["Km2", "v1.V"]\n'
    vnames = m().varnames
    for i,vec in enumerate(dfdp_strings(m, parnames)):
        for j, dxdx in enumerate(vec):
            print '(d d%s/dt / d %s) ='%(vnames[i],parnames[j]), dxdx
    print
    print 'dfdp_strings(m, parnames): (with unknow pars)'
    parnames = "Km3 v1.V".split()
    print 'parnames = ["Km3", "v1.V"]\n'
    vnames = m().varnames
    for i,vec in enumerate(dfdp_strings(m, parnames)):
        for j, dxdx in enumerate(vec):
            print '(d d%s/dt / d %s) ='%(vnames[i],parnames[j]), dxdx
    print
    print '********** Testing _gen_calc_symbmap(m) *******************'
    print '_gen_calc_symbmap(m, with_uncertain = False):'
    print _gen_calc_symbmap(m)
    print '\n_gen_calc_symbmap(m, with_uncertain = True):'
    print _gen_calc_symbmap(m, with_uncertain = True)
    
    print
    print '********** Testing rateCalcString **************************'
    symbmap = _gen_calc_symbmap(m, with_uncertain = False)
    symbmap2 = _gen_calc_symbmap(m, with_uncertain = True)
    for v in (m.v1, m.v2, m.t1, m.t2):
        vstr = v(fully_qualified = True)
        print 'calcstring for %s = %s\n\t'% (get_name(v), vstr), rateCalcString(vstr, symbmap)
    print 'calcstring for v2 with uncertain parameters:\n\t', rateCalcString(m.v2(fully_qualified = True), symbmap2)

    print
    print '********** Testing rate and dXdt generating functions ******'
    print 'Operating point --------------------------------'
    varvalues = [1.0, 0.4]
    pars      = [0.4]
    t         = 0.0
    
    dxdtstrs = [b for (a,b) in dXdt_strings(m)]

    print "t =", t
    print 'variables:'
    pprint.pprint(dict((n, value) for n,value in zip(m().varnames, varvalues)))
    print 'parameters:'
    pprint.pprint(dict((get_name(p), p)     for p in m().parameters))
 
    print '\n---- rates using rates_func(m) -------------------------'
    vratesfunc = rates_func(m)
    vrates = vratesfunc(varvalues,t)
    frmtstr = "%s = %-25s = %s"
    for v,r in zip(m().reactions, vrates):
        print frmtstr % (get_name(v), v(fully_qualified = True), r)

    print '---- transformations using rates_func(m, transf = True) --'
    tratesfunc = rates_func(m,transf = True)
    trates = tratesfunc(varvalues,t)
    for v,r in zip(m().transformations, trates):
        print frmtstr % (get_name(v), v(fully_qualified = True), r)
    print '---- same, at t = 2.0 --'
    trates = tratesfunc(varvalues,2.0)
    for v,r in zip(m().transformations, trates):
        print frmtstr % (get_name(v), v(fully_qualified = True), r)

    print '---- dXdt using getdXdt(m) --------------------------------'
    f = getdXdt(m)
    dXdt = f(varvalues,t)
    for x,s,r in zip(m().varnames, dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x,s,r)

    print '---- dXdt using getdXdt(m) setting uncertain parameters ---'
    print 'f = getdXdt(m, with_uncertain = True)'
    f = getdXdt(m, with_uncertain = True)
    print 'setting uncertain as', dict((get_name(v), value) for v,value in zip(m().uncertain, pars))
    print 'm.set_uncertain(pars)'
    m.set_uncertain(pars)
    dXdt = f(varvalues,t)
    for x,s,r in zip(m().varnames, dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x, s,r)

    print '---- dXdt using dXdt_with(m, pars) ------------------------'
    print 'f = dXdt_with(m, pars)'
    f = dXdt_with(m, pars)
    dXdt   = f(varvalues,t)
    for x,s,r in zip(m().varnames, dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x, s,r)

    print '---- dXdt using getdXdt(m) with a state argument (m.init) --'
    print 'm.init:', m.init
    print 'f = getdXdt(m)'
    f = getdXdt(m)
    print 'dXdt = f(state2array(m,"init"),t)'
    dXdt = f(state2array(m,"init"),t)
    for x,r in zip(m().varnames, dXdt):
        print "d%s/dt = %s" % (x, r)
    print '---- same, changing state argument ---------------------------'
    m.init.A = 2.0
    print 'after m.init.A = 2.0'
    print 'm.init:', m.init
    print 'dXdt = f(state2array(m,"init"),t)'
    dXdt = f(state2array(m,"init"),t)
    for x,r in zip(m().varnames, dXdt):
        print "d%s/dt = %s" % (x, r)
    print '********** Testing add_dSdt_to_model functions ***************'
    print
    m2 = m.clone()
    print m2
    print '----------------------------------------------------'
    pars = "Km2 v1.V".split()
    print 'pars =', pars
    print
    print "!!! applying function add_dSdt_to_model(m, pars) !!!"
    Snames = add_dSdt_to_model(m2, pars)
    print 'Snames = \n', Snames
    print m2
    print '---------------- EXAMPLE 1 ------------------'
    m1 = read_model(models.glyoxalases.text)

    print m1

    solution1 = solve(m1)

    print '--- Last time point ----'
    print 'At t =', solution1.t[-1]
    for x,value in solution1.last:
        print "%-8s= %f" % (x, value)

    #print '---------------- EXAMPLE 2 ------------------'
    m2 = read_model(models.branched.text)

    #print m2
    times = append(linspace(0.0,5.0,500),linspace(5.0,10.0, 500))

    solution2 = solve(m2, tf = 10.0, times=times)

    #print '---------------- EXAMPLE 3 ------------------'
    m3 = read_model(models.ca.text)

    #print m3

    ms = ModelSolver(m3, tf = 8.0, npoints = 2000)
    solution3 = ms.solve()
##     solution3 = solve(m3, tf = 8.0, npoints = 2000)

    print '---------------- EXAMPLE 4 ------------------'
    m4 = read_model(models.rossler.text)

    print m4

    solution4 = solve(m4, tf = 100.0, npoints = 2000, outputs="x1 x2 x3")
    
    def transformation(vars,t):
        if t > 40.0:
            return (vars[0]-5.0, vars[1], vars[2])
        else:
            return (-5.0, vars[1], vars[2])

    solution4.apply_transf(transformation)
    
    #savingfile = open('examples/analysis.png', 'w+b')
    savingfile = 'examples/analysis.png'

    plot ([solution1, solution2, solution3, solution4], superimpose=False, show = True, save2file=savingfile)
    
    if os.path.exists(savingfile):
        print 'removing temp figure file'
        os.remove(savingfile)

if __name__ == "__main__":
    test()
 