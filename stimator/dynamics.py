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
    newlist = [state.varvalues.get(var,0.0) for var in varnames(m)]
    return array(newlist)

def genStoichiometryMatrix(m):
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)

    vnames = varnames(m)
    N = zeros((len(variables(m)),len(reactions(m))), dtype=float)
    for j,v in enumerate(reactions(m)):
        for rORp, signedunit in [(v.reagents,-1.0),(v.products,1.0)]:
            for c in rORp:
                coef, var = (c[1]*signedunit, c[0])
                if var in vnames:
                    ivariable = vnames.index(var) # error handling here
                    N[ivariable, j] = coef
                else:
                    continue # there are no rows for extvariables in stoich. matrix
    return N

def rateCalcString(m, rateString, with_uncertain = False, changing_pars = None):
    # replace varnames
    for i,v in enumerate(variables(m)):
        rateString = re.sub(r"\b"+ v.name+r"\b", "variables[%d]"%i, rateString)

    # replace uncertain parameters or changing parameters
    if with_uncertain:
        for i,u in enumerate(uncertain(m)):
            rateString = re.sub(r"\b"+ u.name+r"\b", "m_Parameters[%d]"%i, rateString)
    else:
        if changing_pars is not None:
            for i,u in enumerate(changing_pars):
                rateString = re.sub(r"\b"+ u+r"\b", "m_Parameters[%d]"%i, rateString)

    # replace parameters
    for p in parameters(m):
        if p.bounds and with_uncertain:
            continue
        rateString = re.sub(r"\b"+ p.name + r"\b", "%g"% p, rateString) 
    return rateString

def rates_strings(m):
    """Generate a tuple of tuples of
       (name, rate) where
       'name' is the name of a reaction
       'rhs' is the string of the rate of the reaction.
    """
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)
    return tuple([(v.name, v.rate) for v in reactions(m)])

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
    for i,x in enumerate(variables(m)):
        name = x.name
        dXdtstring = ''
        for j,v in enumerate(reactions(m)):
            coef = N[i,j]
            if coef == 0.0: continue
            ratestring = '(%s)'% v.rate
            if coef == 1.0:
                ratestring = '+'+ratestring
            else:
                ratestring = "%g*%s" % (coef,ratestring)
                if coef > 0.0:
                    ratestring = '%s%s'%('+', ratestring)
            dXdtstring += ratestring
        res.append((name, dXdtstring))
    return tuple(res)

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
    _dxdtstrings = dXdt_strings(m)
    _symbs = {}
    for x in variables(m):
        _symbs[x.name] = sympy.Symbol(str(x.name))
    for p in parameters(m):
        _symbs[p.name] = sympy.Symbol(str(p.name))
    _nvars = len(variables(m))
    _jfuncs = []
    for _i in range(_nvars):
        _jfuncs.append([])
        _ids = identifiersInExpr(_dxdtstrings[_i][1])
        if len(_ids) == 0:
            for _j in range(_nvars):
                _jfuncs[_i].append('0.0')
        else:
        
            for _j in range(_nvars):
                _res = eval(_dxdtstrings[_i][1], None, _symbs)
                _res = _res * _scale
                _dres = str(sympy.diff(_res, _symbs[variables(m)[_j].name]))
                if _dres == '0':
                    _dres = '0.0'
                _jfuncs[_i].append(_dres)
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
        print 'ERROR: sympy module must be installed to generate partial derivative strings'
        raise
    _dxdtstrings = dXdt_strings(m)
    _symbs = {}
    for x in variables(m):
        _symbs[x.name] = sympy.Symbol(str(x.name))
    for p in parameters(m):
        _symbs[p.name] = sympy.Symbol(str(p.name))
    _nvars = len(variables(m))
    _npars = len(_parnames)
    _jfuncs = []
    for _i in range(_nvars):
        _jfuncs.append([])
        _ids = identifiersInExpr(_dxdtstrings[_i][1])
        if len(_ids) == 0:
            for _j in range(_npars):
                _jfuncs[_i].append('0.0')
        else:
        
            for _j in range(_npars):
                _res = eval(_dxdtstrings[_i][1], None, _symbs)
                _res = _res * _scale
                if not _symbs.has_key(_parnames[_j]):
                    _dres = '0.0'
                else:
                    _dres = str(sympy.diff(_res, _symbs[_parnames[_j]]))
                if _dres == '0':
                    _dres = '0.0'
                _jfuncs[_i].append(_dres)
                #~ print _jfuncs
    return _jfuncs
        
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
    ratebytecode = [compile(rateCalcString(m, v.rate, with_uncertain=with_uncertain), '<string>','eval') for v in collection]
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
    for a in argnames:
        i, collection = m._Model__findComponent(a)
        if collection == 'parameters':
            data.append(('p', getattr(self, a)))
        elif collection == 'variables':
            data.append(('v', i))
        elif collection == 'transf':
            data.append(('t', compile(rateCalcString(m, transformations(m)[i].rate), '<string>','eval')))
        elif collection == 'reactions':
            data.append(('r', compile(rateCalcString(m, reactions(m)[i].rate), '<string>','eval')))
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

## def getdXdt(m, with_uncertain = False, scale = 1.0, t0=0.0, changing_pars = None):
##     """Generate function to compute rhs of SODE for this model.
##     
##        Function has signature f(variables, t)
##        This is compatible with scipy.integrate.odeint"""

##     check, msg = m.checkRates()
##     if not check:
##         print "vars = "
##         print [x.name for x in variables(m)]
##         raise BadRateError(msg)
##     #compile rhs
##     
##     rhsides = dXdt_strings(m)
##     ratebytecode = [compile(rateCalcString(m, "%g *(%s)"%(scale,rhs), 
##                                            with_uncertain = with_uncertain, 
##                                            changing_pars=changing_pars), 
##                                            '<string>','eval') for (xname,rhs) in rhsides]
##     # create array to hold v's
##     x = empty(len(variables(m)))
##     en = list(enumerate(ratebytecode))
##     
##     def f2(variables, t):
##         m_Parameters = m._Model__m_Parameters
##         t = t*scale + t0
##         for i,r in en:
##             x[i] = eval(r)
##         return x
##     return f2

def getdXdt(m, with_uncertain = False, scale = 1.0, t0=0.0, changing_pars = None):
    """Generate function to compute rhs of SODE for this model.
    
       Function has signature f(variables, t)
       This is compatible with scipy.integrate.odeint"""

    check, msg = m.checkRates()
    if not check:
        print "vars = "
        print [x.name for x in variables(m)]
        raise BadRateError(msg)
    #compile rate laws
    ratebytecode = [compile(rateCalcString(m, v.rate, 
                                           with_uncertain = with_uncertain, 
                                           changing_pars=changing_pars), 
                                           '<string>','eval') for v in reactions(m)]
    # compute stoichiometry matrix, scale and transpose
    N  = genStoichiometryMatrix(m)
    N *= scale
    NT = N.transpose()
    # create array to hold v's
    v = empty(len(reactions(m)))
    x = empty(len(variables(m)))
    en = list(enumerate(ratebytecode))
    
    def f2(variables, t):
        m_Parameters = m._Model__m_Parameters
        t = t*scale + t0
        for i,r in en:
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
    ratebytecode = [compile(rateCalcString(m, v.rate, with_uncertain = True), '<string>','eval') for v in reactions(m)]
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
    ratebytecode = [[compile(rateCalcString(m, col, with_uncertain = with_uncertain), '<string>','eval') for col in line] for line in Jstrings]

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
    v1 = A -> B, rate = V1*A/(Km1 + A)
    v2 = B ->  , rate = V2*B/(Km2 + B)
    V1  = 1
    Km1 = 1
    V2  = sqrt(4.0)
    Km2 = 0.2
    find Km2 in [0, 1.2]
    init = state(B = 0.4, A = 1)
    ~ t1 = A+B
    ~ t2 = V1*A * step(t, 1.0)
    """)
    print m

    print '********** Testing rate and dXdt strings *******************'
    print 'rates_strings(): -------------------------'
    for v in rates_strings(m):
        print v
    print '\ndXdt_strings(): -------------------------'
    for dxdt in dXdt_strings(m):
        print dxdt
    print
    print '********** Testing stoichiometry matrix ********************'
    print 'Stoichiometry matrix:'
    N = genStoichiometryMatrix(m)
    print '  ', '  '.join([v.name for v in reactions(m)])
    for i,x in enumerate(variables(m)):
        print x.name, N[i, :]
    print
    print '********** Testing state2array()****************************'
    print 'vectorize(m,"init"):'
    v = state2array(m,"init")
    print v, 'of type', type(v)
    print
    print '********** Testing rateCalcString **************************'
    print 'calcstring for v1:\n', rateCalcString(m, m.v1.rate)
    print
    print 'calcstring for v2:\n', rateCalcString(m, m.v2.rate)
    print
    print 'calcstring for v2 with uncertain parameters:\n', rateCalcString(m, m.v2.rate, True)
    print
    print 'calcstring for t1:\n', rateCalcString(m, m.t1.rate)
    print
    print 'calcstring for t2:\n', rateCalcString(m, m.t2.rate)
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
    pprint.pprint(dict((p.name, p)     for p in parameters(m)))
 
    print '---- rates using rates_func(m) -------------------------'
    vratesfunc = rates_func(m)
    vrates = vratesfunc(varvalues,t)
    for v,r in zip(reactions(m), vrates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)

    print '---- transformations using rates_func(m, transf = True) --'
    tratesfunc = rates_func(m,transf = True)
    trates = tratesfunc(varvalues,t)
    for v,r in zip(transformations(m), trates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)
    print '---- same, at t = 2.0 --'
    trates = tratesfunc(varvalues,2.0)
    for v,r in zip(transformations(m), trates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)

    print '---- dXdt using getdXdt(m) --------------------------------'
    f = getdXdt(m)
    dXdt = f(varvalues,t)
    for x,s,r in zip(variables(m), dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x.name, s,r)

    print '---- dXdt using getdXdt(m) setting uncertain parameters ---'
    print 'f = getdXdt(m, with_uncertain = True)'
    f = getdXdt(m, with_uncertain = True)
    print 'setting uncertain as', dict((v.name, value) for v,value in zip(uncertain(m), pars))
    print 'm.set_uncertain(pars)'
    m.set_uncertain(pars)
    dXdt = f(varvalues,t)
    for x,s,r in zip(variables(m), dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x.name, s,r)

    print '---- dXdt using dXdt_with(m, pars) ------------------------'
    print 'f = dXdt_with(m, pars)'
    f = dXdt_with(m, pars)
    dXdt   = f(varvalues,t)
    for x,s,r in zip(variables(m), dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x.name, s,r)

    print '---- dXdt using getdXdt(m) with a state argument (m.init) --'
    print 'm.init:', m.init
    print 'f = getdXdt(m)'
    f = getdXdt(m)
    print 'dXdt = f(state2array(m,"init"),t)'
    dXdt = f(state2array(m,"init"),t)
    for x,r in zip(variables(m), dXdt):
        print "d%s/dt = %s" % (x.name, r)
    print '---- same, changing state argument ---------------------------'
    m.init.A = 2.0
    print 'after m.init.A = 2.0'
    print 'm.init:', m.init
    print 'dXdt = f(state2array(m,"init"),t)'
    dXdt = f(state2array(m,"init"),t)
    for x,r in zip(variables(m), dXdt):
        print "d%s/dt = %s" % (x.name, r)

if __name__ == "__main__":
    test()
 