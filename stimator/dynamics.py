#!/usr/bin/env python
# -*- coding: latin1-*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator dynamical systems related functions
# Copyright António Ferreira 2006-2011
#----------------------------------------------------------------------------
from model import *
import math
from kinetics import *
from numpy import *
import pprint
from modelparser import read_model

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
    return array([state._ownparameters.get(var,0.0) for var in varnames(m)])

def genStoichiometryMatrix(m):
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)

    vnames = varnames(m)
    N = zeros((len(vnames),len(reactions(m))), dtype=float)
    for j,v in enumerate(reactions(m)):
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
    return tuple([(get_name(v), v(fully_qualified = fully_qualified)) for v in reactions(m)])

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
    for i,name in enumerate(varnames(m)):
        dXdtstring = ''
        for j,v in enumerate(reactions(m)):
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
    for x in varnames(m):
        symbname = '_symbol_Id%d'% symbcounter
        symbmap[x] = symbname
        sympysymbs[symbname] = sympy.Symbol(symbname)
        symbcounter += 1
    for p in parameters(m):
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
    _dxdtstrings = [_d[1] for _d in dXdt_strings(m)]
    _nvars = len(_dxdtstrings)
    _varnames = varnames(m)
    _symbmap, _sympysymbs = _gen_canonical_symbmap(m)
    for _i in range(_nvars):
        _dxdtstrings[_i] = _replace_exprs2canonical(_dxdtstrings[_i], _symbmap)

    _jfuncs = []
    for _i in range(_nvars):
        _jfuncs.append([])
        _ids = identifiersInExpr(_dxdtstrings[_i])
        if len(_ids) == 0:
            for _j in range(_nvars):
                _jfuncs[_i].append('0.0')
        else:
            for _j in range(_nvars):
                _varsymb = _symbmap[_varnames[_j]]
                _res = eval(_dxdtstrings[_i], None, _sympysymbs)
                _res = _res * _scale
                _dres = str(sympy.diff(_res, _varsymb))
                if _dres == '0':
                    _dres = '0.0'
                _jfuncs[_i].append(_dres)
    # back to original ids
    for _i in range(_nvars):
        for _j in range(_nvars):
            _jfuncs[_i][_j] = _replace_canonical2exprs(_jfuncs[_i][_j], _symbmap)
    return _jfuncs
        
def dfdp_strings(m, _parnames, _scale = 1.0):
    """Generate a matrix (list of lists) of strings
       to compute the partial derivatives of rhs of SODE
       with respect to a list of parameters.
       _parnames is a list of parameter names.
    
       IMPORTANT: sympy module must be installed!"""

    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    try:
        import sympy
    except:
        print 'ERROR: sympy module must be installed to generate Jacobian strings'
        raise
    _dxdtstrings = [_d[1] for _d in dXdt_strings(m)]
    _nvars = len(_dxdtstrings)
    _symbmap, _sympysymbs = _gen_canonical_symbmap(m)
    for _i in range(_nvars):
        _dxdtstrings[_i] = _replace_exprs2canonical(_dxdtstrings[_i], _symbmap)

    _npars = len(_parnames)
    _jfuncs = []
    for _i in range(_nvars):
        _jfuncs.append([])
        _ids = identifiersInExpr(_dxdtstrings[_i])
        if len(_ids) == 0:
            for _j in range(_npars):
                _jfuncs[_i].append('0.0')
        else:
            for _j in range(_npars):
                if _parnames[_j] not in _symbmap:
                    _dres = '0.0'
                else:
                    _varsymb = _symbmap[_parnames[_j]]
                    _res = eval(_dxdtstrings[_i], None, _sympysymbs)
                    _res = _res * _scale
                    _dres = str(sympy.diff(_res, _varsymb))
                    if _dres == '0':
                        _dres = '0.0'
                _jfuncs[_i].append(_dres)
    # back to original ids
    for _i in range(_nvars):
        for _j in range(_npars):
            _jfuncs[_i][_j] = _replace_canonical2exprs(_jfuncs[_i][_j], _symbmap)
    return _jfuncs
        
def _gen_calc_symbmap(m, with_uncertain = False):
    symbmap = {}
    for i, x in enumerate(varnames(m)):
        symbname = "variables[%d]"%i
        symbmap[x] = symbname
    if with_uncertain:
        for i,u in enumerate(uncertain(m)):
            symbname = "m_Parameters[%d]"%i
            symbmap[get_name(u)] = symbname
    for p in parameters(m):
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
        collection = transformations(m)
    else:
        collection = reactions(m)

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
        if not isinstance(f,str):
            raise TypeError('argument must be a string or a callable.')
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
            vstr = rateCalcString(transformations(m)[i](fully_qualified = True), symbmap)
            data.append(('t', compile(vstr, '<string>','eval')))
        elif collection == 'reactions':
            vstr = rateCalcString(reactions(m)[i](fully_qualified = True), symbmap)
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
    ratebytecode = compile_rates(m, reactions(m), with_uncertain = with_uncertain)

    # compute stoichiometry matrix, scale and transpose
    N  = genStoichiometryMatrix(m)
    N *= scale
    NT = N.transpose()
    # create array to hold v's
    v = empty(len(reactions(m)))
    x = empty(len(varnames(m)))
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
    ratebytecode = compile_rates(m, reactions(m), with_uncertain = with_uncertain)
    # compute stoichiometry matrix, scale and transpose
    N  = genStoichiometryMatrix(m)
    N *= scale
    NT = N.transpose()
    # create array to hold v's
    v = empty(len(reactions(m)))
    x = empty(len(varnames(m)))
    for i in range(len(reactions(m))):
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
    ratebytecode = compile_rates(m, reactions(m), with_uncertain = True)
    # compute stoichiometry matrix, scale and transpose
    N  = genStoichiometryMatrix(m)
    N *= scale
    NT = N.transpose()

    # create array to hold v's
    v = empty(len(reactions(m)))
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

def test():
     
    m = read_model("""
    title a simple 2 enzyme system
    v1 = A -> B, rate = V*A/(Km + A), V = 1, Km = 1
    v2 = B ->  , rate = V*B/(Km2 + B)
    V  = sqrt(4.0)
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
    print '  ', '  '.join([get_name(v) for v in reactions(m)])
    for i,x in enumerate(varnames(m)):
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
    vnames = varnames(m)
    for i,vec in enumerate(Jacobian_strings(m)):
        for j, dxdx in enumerate(vec):
            print '(d d%s/dt / d %s) ='%(vnames[i],vnames[j]), dxdx
    print
    print 'dfdp_strings(m, parnames): ------------------'
    parnames = "Km2 v1.V".split()
    print '\nparnames = ["Km2", "v1.V"]\n'
    vnames = varnames(m)
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
    pprint.pprint(dict((n, value) for n,value in zip(varnames(m), varvalues)))
    print 'parameters:'
    pprint.pprint(dict((get_name(p), p)     for p in parameters(m)))
 
    print '---- rates using rates_func(m) -------------------------'
    vratesfunc = rates_func(m)
    vrates = vratesfunc(varvalues,t)
    frmtstr = "%s = %-25s = %s"
    for v,r in zip(reactions(m), vrates):
        print frmtstr % (get_name(v), v(fully_qualified = True), r)

    print '---- transformations using rates_func(m, transf = True) --'
    tratesfunc = rates_func(m,transf = True)
    trates = tratesfunc(varvalues,t)
    for v,r in zip(transformations(m), trates):
        print frmtstr % (get_name(v), v(fully_qualified = True), r)
    print '---- same, at t = 2.0 --'
    trates = tratesfunc(varvalues,2.0)
    for v,r in zip(transformations(m), trates):
        print frmtstr % (get_name(v), v(fully_qualified = True), r)

    print '---- dXdt using getdXdt(m) --------------------------------'
    f = getdXdt(m)
    dXdt = f(varvalues,t)
    for x,s,r in zip(varnames(m), dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x,s,r)

    print '---- dXdt using getdXdt(m) setting uncertain parameters ---'
    print 'f = getdXdt(m, with_uncertain = True)'
    f = getdXdt(m, with_uncertain = True)
    print 'setting uncertain as', dict((get_name(v), value) for v,value in zip(uncertain(m), pars))
    print 'm.set_uncertain(pars)'
    m.set_uncertain(pars)
    dXdt = f(varvalues,t)
    for x,s,r in zip(varnames(m), dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x, s,r)

    print '---- dXdt using dXdt_with(m, pars) ------------------------'
    print 'f = dXdt_with(m, pars)'
    f = dXdt_with(m, pars)
    dXdt   = f(varvalues,t)
    for x,s,r in zip(varnames(m), dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x, s,r)

    print '---- dXdt using getdXdt(m) with a state argument (m.init) --'
    print 'm.init:', m.init
    print 'f = getdXdt(m)'
    f = getdXdt(m)
    print 'dXdt = f(state2array(m,"init"),t)'
    dXdt = f(state2array(m,"init"),t)
    for x,r in zip(varnames(m), dXdt):
        print "d%s/dt = %s" % (x, r)
    print '---- same, changing state argument ---------------------------'
    m.init.A = 2.0
    print 'after m.init.A = 2.0'
    print 'm.init:', m.init
    print 'dXdt = f(state2array(m,"init"),t)'
    dXdt = f(state2array(m,"init"),t)
    for x,r in zip(varnames(m), dXdt):
        print "d%s/dt = %s" % (x, r)

if __name__ == "__main__":
    test()
 