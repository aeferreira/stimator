import re
import math
from itertools import chain

import numpy as np
from scipy import integrate
import sympy

from stimator.timecourse import SolutionTimeCourse, Solutions
from stimator.utils import _is_string, _is_sequence

from stimator.examples import models

#-----------------------------------------
# BadRateError Exception
#-----------------------------------------

class BadRateError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)

def identifiersInExpr(_expr):
    iterator = identifier.finditer(_expr)
    return [_expr[m.span()[0]:m.span()[1]] for m in iterator]


def init2array(model):
    """Transforms a state object into a numpy.array object.

       This is necessary for most numerical functions of numpy+scipy.
       Can accept the name of a state (must exist in Model) or state object.
       Values are returned in the order of model variables.
    """
    return np.array([model.get_init(var) for var in model.varnames])


def genStoichiometryMatrix(model):
    check, msg = model.checkRates()
    if not check:
        raise BadRateError(msg)

    vnames = model.varnames
    N = np.zeros((len(vnames), len(model.reactions)), dtype=float)
    for j, v in enumerate(model.reactions):
        for rORp, sign_one in [(v._reagents, -1.0), (v._products, 1.0)]:
            for var, coef in rORp:
                if var in vnames:
                    ivar = vnames.index(var)
                    N[ivar, j] = coef * sign_one
                else:
                    continue # no rows for extvariables.
    return N


def rates_strings(model, fully_qualified=True):
    """Generate a dictionary of name: rate.

       'name' is the name of a reaction
       'rate' is the string of the rate of the reaction.
    """
    check, msg = model.checkRates()
    if not check:
        raise BadRateError(msg)
    return {v.name: v(fully_qualified=fully_qualified) for v in model.reactions}

def dXdt_strings(model):
    """Generate a dictionary of name: rhs.
       
       'name' is the name of a variable
       'rhs' is the string of the rhs of that variable in the SODE.
    """

    check, msg = model.checkRates()
    if not check:
        raise BadRateError(msg)
    symbols = _gen_canonical_symbmap(model)

    N = genStoichiometryMatrix(model)
    res = {}
    for i, name in enumerate(model.varnames):
        dXdtstring = ''
        for j, v in enumerate(model.reactions):
            coef = N[i, j]
            if coef == 0.0:
                continue
            ratestring = '(%s)'% v(fully_qualified=True)
            if coef == 1.0:
                ratestring = '+'+ ratestring
            else:
                ratestring = "%g*%s" % (coef, ratestring)
                if coef > 0.0:
                    ratestring = '%s%s'%('+', ratestring)
            dXdtstring += ratestring
        dXdtstring = _simplify_expr(dXdtstring, symbols)
        res[name] = dXdtstring
    return res


def _gen_canonical_symbmap(m, extra_id_list=None):
    check, msg = m.checkRates()
    if not check:
        raise BadRateError(msg)

    symbmap = {}
    sympysymbs = {}

    par_names = [p.name for p in m.parameters]
    if extra_id_list is None:
        extra_id_list = []
    
    for i, x in enumerate(chain(m.varnames, par_names, extra_id_list)):
        name = '_symbol_Id%d'% i
        symbmap[x] = name
        sympysymbs[name] = sympy.Symbol(name)

    return {'s_table': symbmap, 'sympy_s_table': sympysymbs}

def _replace_exprs2canonical(s, symbmap):
    for symb in symbmap:
        symbesc = symb.replace('.', '\.')
        s = re.sub(r"(?<![_.])\b%s\b(?![_.\[])"%symbesc, symbmap[symb], s)
    return s

def _replace_canonical2exprs(s, symbmap):
    for symb in symbmap:
        s = re.sub(r"(?<![.])\b%s\b" % symbmap[symb], symb, s)
    return s

def _differentiate_expr(expr, wrt, symbols, _scale=1.0):
    symbmap, sympysymbs = symbols['s_table'], symbols['sympy_s_table']

    texpr = _replace_exprs2canonical(expr, symbmap)
    if wrt not in symbmap:
        return '0.0'
    varsymb = symbmap[wrt]
    ids = identifiersInExpr(texpr)
    if len(ids) == 0 or varsymb not in ids:
        return '0.0'

    res = eval(texpr, None, sympysymbs)
    if _scale != 1.0:
        res = res * _scale
    dres = str(sympy.diff(res, varsymb))
    if dres == '0':
        dres = '0.0'
    dres = _replace_canonical2exprs(dres, symbmap)
    return dres

def _simplify_expr(expr, symbols):
    symbmap, sympysymbs = symbols['s_table'], symbols['sympy_s_table']
    resstr = _replace_exprs2canonical(expr, symbmap)
    #make sympy reduce the expression using sympysymbs dictionary
    dres = str(eval(resstr, None, sympysymbs))
    if dres == '0':
        dres = '0.0'
    dres = _replace_canonical2exprs(dres, symbmap)
    return dres


def Jacobian_strings(m, _scale=1.0, symbols=None):
    """Generate a matrix (list of lists) of strings
       to compute the jacobian for this model."""

    dxdtstrings = dXdt_strings(m)

    x_names = m.varnames

    if symbols is None:
        symbols = _gen_canonical_symbmap(m)

    jfuncs = []
    for x in x_names:
        dlist = []
        for y in x_names:
            dexpr = _differentiate_expr(dxdtstrings[x], y, symbols, _scale)
            dlist.append(dexpr)
        jfuncs.append(dlist)
    return jfuncs


def dfdp_strings(m, parnames, _scale=1.0, symbols=None):
    """Generate a matrix (list of lists) of strings
       to compute the partial derivatives of rhs of SODE
       with respect to a list of parameters.
       parnames is a list of parameter names."""

    dxdtstrings = dXdt_strings(m)

    vnames = m.varnames

    if symbols is None:
        symbols = _gen_canonical_symbmap(m)

    dxdp_strs = []
    for x in vnames:
        dlist = []
        for p in parnames:
            dexpr = _differentiate_expr(dxdtstrings[x], p,
                                        symbols, _scale)
            dlist.append(dexpr)
        dxdp_strs.append(dlist)
    return dxdp_strs


def add_dSdt_to_model(m, pars):
    """Add sensitivity ODEs to model, according to formula:

    dS/dt = df/dx * S + df/dp

    m is a model object
    pars are a list of parameter names
    """

    #Find pars that are initial values
    init_of = []
    for p in pars:
        init_of.append(None)
        if '.' in p:
            tks = p.split('.')
            if tks[1] in m.init:
                init_of[-1] = tks[1]

    #create names and symbols for sensitivities
    Snames = []
    Smatrix = []
    vnames = m.varnames
    for x in vnames:
        new_row = []
        for p in pars:
            Sname = "d_%s_d_%s" % (x, p.replace('.', '_'))
            new_row.append(Sname)
            Snames.append((x, p, Sname))
        Smatrix.append(new_row)

    id_list = [name[2] for name in Snames]

    symbols = _gen_canonical_symbmap(m, extra_id_list=id_list)

    J = Jacobian_strings(m, symbols=symbols)
    dfdpstrs = dfdp_strings(m, pars, symbols=symbols)
    nvars = len(J)
    npars = len(pars)

    #compute rhs of sensitivities symbolically
    for i in range(nvars):
        x = m.varnames[i]
        for j in range(npars):
            #compute string for dS/dt
            resstr = dfdpstrs[i][j] if init_of[j] is None else ''
            # matrix multiplication with strings:
            for k in range(nvars):
                resstr = resstr + "+(%s)*(%s)"%(J[i][k], Smatrix[k][j])
            dres = _simplify_expr(resstr, symbols=symbols)

            # add new variable to model
            name = Smatrix[i][j]
            m.set_variable_dXdt(name, dres)

            # set initial value
            if init_of[j] is None:
                m.set_init([(name, 0.0)])
            else:
                if init_of[j] == x:
                    m.set_init([(name, 1.0)])
                else:
                    m.set_init([(name, 0.0)])
    return Snames

def rateCalcString(rateString, symbmap):
    return _replace_exprs2canonical(rateString, symbmap)

def _gen_calc_symbmap(model, with_uncertain=False):
    symbmap = {}
    for i, x in enumerate(model.varnames):
        symbname = "variables[%d]"%i
        symbmap[x] = symbname
    for i, invar in enumerate(model.input_variables):
        symbname = "input_variables[%d]"%i
        symbmap[invar.name] = symbname
    if with_uncertain:
        for i, u in enumerate(model.with_bounds):
            symbname = "m_Parameters[%d]"%i
            symbmap[u.name] = symbname
    for p in model.parameters:
        if p.bounds and with_uncertain:
            continue
        valuestr = "%g" % p
        symbmap[p.name] = valuestr

    return symbmap


def compile_dxdt(model, with_uncertain=False):
    dxdtstrings = dXdt_strings(model)
    symbmap = _gen_calc_symbmap(model, with_uncertain=with_uncertain)
    dxdt_exprs = [rateCalcString(expr, symbmap) for expr in dxdtstrings]
    return [compile(expr, '<string>', 'eval') for expr in dxdt_exprs]


def compile_all_rates(model, with_uncertain=False):
    symbmap = _gen_calc_symbmap(model, with_uncertain=with_uncertain)
    input_rates = [rateCalcString(v(fully_qualified=True), symbmap) for v in model.input_variables]
    rate_rates = [rateCalcString(v(fully_qualified=True), symbmap) for v in model.reactions]
    transf_rates = [rateCalcString(v(fully_qualified=True), symbmap) for v in  model.transformations]
    input_bc = [compile(v, '<string>', 'eval') for v in input_rates]
    rate_bc = [compile(v, '<string>', 'eval') for v in rate_rates]
    transf_bc = [compile(v, '<string>', 'eval') for v in transf_rates]
    return input_bc, rate_bc, transf_bc


def _get_rates_function(model, with_uncertain):
    check, msg = model.checkRates()
    if not check:
        raise BadRateError(msg)

    # compile all changing variables
    input_bc, rate_bc, transf_bc = compile_all_rates(model, with_uncertain=with_uncertain)

    # create arrays to hold computed values
    input_variables = np.empty(len(model.input_variables))
    _v_rates = np.empty(len(model.reactions))
    _t_rates = np.empty(len(model.transformations))

    def f(variables, t):
        m_Parameters = model._Model__m_Parameters
##         print('************** locals *********************')
##         for l, v in locals().items():
##             print(l,'--->', v)
##         print('************** end locals *********************')

        for i, r in enumerate(input_bc):
            input_variables[i] = eval(r, model._usable_functions, locals())
        for i, r in enumerate(rate_bc):
            _v_rates[i] = eval(r, model._usable_functions, locals())
        for i, r in enumerate(transf_bc):
            _t_rates[i] = eval(r, model._usable_functions, locals())
        return input_variables, _v_rates, _t_rates

    return f


def get_outputs_function(model, with_uncertain=False, out_names=None):
    check, msg = model.checkRates()
    if not check:
        raise BadRateError(msg)

    # compile all changing variables
    input_bc, rate_bc, transf_bc = compile_all_rates(model, with_uncertain=with_uncertain)

    # create arrays to hold computed values
    input_variables = np.empty(len(model.input_variables))

    # handle names of output variables to retain
    #loc_codes = [_get_code_index(model, name) for name in out_names]

    out_bytecode = []

    for name in out_names:
        exist = False

        for i, p in enumerate(model.with_bounds):
            if p.name == name:
                out_bytecode.append(compile('m_Parameters[{}]'.format(i), '<string>', 'eval'))
                exist = True
                break
        if exist:
            continue
        for p in model.parameters:
            if p.name == name:
                out_bytecode.append(compile(str(model.getp(name)), '<string>', 'eval'))
                exist = True
                break
        if exist:
            continue
        i = model._Model__reactions.iget(name)
        if i is not None:
            out_bytecode.append(rate_bc[i])
            continue
        try:
            i = model._Model__variables.index(name)
            out_bytecode.append(compile('variables[{}]'.format(i), '<string>', 'eval'))
            continue
        except:
            pass
        i = model._Model__transf.iget(name)
        if i is not None:
            out_bytecode.append(transf_bc[i])
            continue
        i = model._Model__invars.iget(name)
        if i is not None:
            out_bytecode.append(compile('input_variables[{}]'.format(i), '<string>', 'eval'))
            continue
        if not exist:
            raise AttributeError('{} is not a component in this model'.format(name))

    _out_rates = np.empty(len(out_bytecode))

    def out_f(variables, t):
        m_Parameters = model._Model__m_Parameters

        for i, r in enumerate(input_bc):
            input_variables[i] = eval(r, model._usable_functions, locals())
        for i, r in enumerate(out_bytecode):
            _out_rates[i] = eval(r, model._usable_functions, locals())
        return _out_rates

    return out_f


def all_rates_func(m, with_uncertain=False, scale=1.0, t0=0.0):
    """Generate function to compute rate vector for this model.

       Function has signature f(variables, t)"""

    get_rates = _get_rates_function(m, with_uncertain=with_uncertain)

    def fout(variables, t):
        t = t*scale + t0
        return get_rates(variables, t)

    return fout


def getdXdt(model, with_uncertain=False, scale=1.0, t0=0.0):
    """Generate function to compute rhs of SODE for this model.

       Function has signature f(variables, t)
       This is compatible with scipy.integrate.odeint"""

    f_rates = _get_rates_function(model, with_uncertain=with_uncertain)

    # compute stoichiometry matrix, scale and transpose
    N = genStoichiometryMatrix(model)
    N *= scale
    NT = N.transpose()
    dxdt = np.empty(len(model.varnames))

    def fout(variables, t):
        t = t*scale + t0
        _, v, _ = f_rates(variables, t)
        np.dot(v, NT, dxdt)
        return dxdt
    return fout


def getJacobian(m, with_uncertain=False, scale=1.0, t0=0.0):
    """Generate function to compute the jacobian for this model.

       Function has signature J(variables, t)
       and returns an nvars x nvars numpy array
       IMPORTANT: sympy module must be installed!"""

    Jstrings = Jacobian_strings(m, _scale=scale)
    nvars = len(Jstrings)

    #compile rate laws
    symbmap = _gen_calc_symbmap(m, with_uncertain=with_uncertain)
    ratestrs = [[rateCalcString(col, symbmap) for col in line] for line in Jstrings]
    ratebytecode = [[compile(col, '<string>', 'eval') for col in line] for line in ratestrs]

    def J(variables, t):
        Jarray = np.empty((nvars, nvars), float)
        t = t*scale + t0
        for i in range(nvars):
            for j in range(nvars):
                Jarray[i, j] = eval(ratebytecode[i][j], m._usable_functions, locals())
        return Jarray
    def J2(variables, t):
        m_Parameters = m._Model__m_Parameters
        Jarray = np.empty((nvars, nvars), float)
        t = t*scale + t0
        for i in range(nvars):
            for j in range(nvars):
                Jarray[i, j] = eval(ratebytecode[i][j], m._usable_functions, locals())
        return Jarray
    if with_uncertain:
        return J2
    else:
        return J


def get_outputs_decl(model, ignore_replist=False):
    decl = model.metadata.get('!!', None)

    if decl is not None and not ignore_replist:
        names = decl.strip().split()
        return [n.strip() for n in names]
    else:
        return []


def process_outputs_list(names, model):
    if names is None:
        return None
    if _is_string(names):
        names = [names.strip()]
    special_transf = ['~']
    special_rates = ['>', '>>', '->']

    if not _is_sequence(names):
        raise TypeError('outputs must be a sequence of names.')

    if len(names) == 0:
        return None

    out_names = []
    for a in names:
        if not _is_string(a):
            raise TypeError(str(a) + ' must be a string')
        if a in special_transf:
            out_names.extend([x.name for x in model.transformations])
        elif a in special_rates:
            out_names.extend([x.name for x in model.reactions])
        else:
            out_names.append(a)
    return out_names


def solve(model,
          tf=None,
          npoints=500,
          t0=0.0,
          initial=None,
          times=None,
          outputs=None,
          title=None,
          ignore_replist=False):

    #solver = integrate._odepack.odeint
    solver = integrate.odeint
    names = [x for x in model.varnames]

    #get initial values
    if initial is None:
        y0 = init2array(model)
    else:
        y0 = np.copy(initial)
    if tf is None:
        try:
            tf = float(model.metadata.get('tf', None))
        except:
            pass
        if tf is None:
            tf = 1.0
    if times is None:
        times = np.linspace(t0, tf, npoints)

    # scale times to maximum time in data
    t0 = times[0]
    scale = float(times[-1] - t0)
    #scale = 1.0

    f = getdXdt(model, scale=scale, t0=t0)
    transf_f = None

    # get outputs
    # overide if outputs argument is not None
    out_names = get_outputs_decl(model, ignore_replist)
    if outputs is not None:
        out_names = outputs
    out_names = process_outputs_list(out_names, model)

    if out_names is not None:
        transf_f = get_outputs_function(model, out_names=out_names)
        transf_f.names = out_names

    t = np.copy((times-t0)/scale)  # this scales time points

    output = solver(f, y0, t,
                    args=(),
                    Dfun=None,
                    col_deriv=0,
                    full_output=True,
                    ml=None,
                    rtol=None,
                    mu=None,
                    atol=None,
                    tcrit=None,
                    h0=0.0,
                    hmax=0.0,
                    hmin=0.0,
                    ixpr=0,
                    mxstep=0,
                    mxhnil=0,
                    mxordn=12,
                    mxords=5)#, tfirst=False)
    out_message = output[1]['message'].strip()
    if out_message != 'Integration successful.':
        print('Solution failed:', out_message)
        return None

    Y = output[0]

    if title is None:
        title = model.metadata.get('title', '')
    Y = np.copy(Y.T)

    sol = SolutionTimeCourse(times, Y, names, title, dense=True)

    if transf_f is not None:
        sol.apply_transf(transf_f, transf_f.names)

    return sol

class ModelSolver(object):
    def __init__(self,
                 model,
                 tf=1.0,
                 npoints=500,
                 t0=0.0,
                 initial=None,
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
        if initial is None:
            self.y0 = np.copy(init2array(self.model))
        else:
            self.y0 = np.copy(initial)

        self.times = times
        if self.times is None:
            self.times = np.linspace(t0, tf, npoints)

        # scale times to maximum time in data
        t0 = self.times[0]
        scale = float(self.times[-1] - t0)
        self.t = (self.times-t0)/scale  # this scales time points

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
                mapinit2pars.append((ix, i))
            self.model.set_bounds(parname, (0, 1)) # bogus bounds

        self.pars_initindexes = np.array([j for (i, j) in mapinit2pars], dtype=int)
        self.vars_initindexes = np.array([i for (i, j) in mapinit2pars], dtype=int)

        self.f = getdXdt(self.model, with_uncertain=True, scale=scale, t0=t0)

        self.transf_f = None
        # get outputs
        # overide if outputs argument is not None
        out_names = get_outputs_decl(self.model, ignore_replist)
        if outputs is not None:
            out_names = outputs
        out_names = process_outputs_list(out_names, model)

        if out_names is not None:
            self.transf_f = get_outputs_function(model, with_uncertain=True, out_names=out_names)
            self.transf_f.names = out_names

    def solutions_names(self):
        if self.transf_f is not None:
            return self.transf_f.names
        return self.model.varnames

    def solve(self, title=None, par_values=None, npoints=None):

        # set initial values
        y0 = np.copy(self.y0)

        # set varying parameters (may be initial values)
        if par_values is not None:
            par_values = np.array(par_values)
            self.model.set_uncertain(par_values)
            # fill uncertain initial values
            y0[self.vars_initindexes] = par_values[self.pars_initindexes]

        if npoints is not None:
            tpoints = np.linspace(self.t[0], self.t[-1], npoints)
        else:
            tpoints = self.t


        output = integrate.odeint(self.f, y0, tpoints,
                                  args=(),
                                  Dfun=None,
                                  col_deriv=0,
                                  full_output=True,
                                  ml=None,
                                  rtol=None,
                                  mu=None,
                                  atol=None,
                                  tcrit=None,
                                  h0=0.0,
                                  hmax=0.0,
                                  hmin=0.0,
                                  ixpr=0,
                                  mxstep=0,
                                  mxhnil=0,
                                  mxordn=12,
                                  mxords=5)#, tfirst=False)
        out_message = output[1]['message'].strip()
        if out_message != 'Integration successful.':
            print('Solution failed:', out_message)
            return None

        Y = output[0]

        if title is None:
            title = self.title

        if npoints is not None:
            tpoints = np.linspace(self.times[0], self.times[-1], npoints)
        else:
            tpoints = self.times

        sol = SolutionTimeCourse(tpoints, Y.T, self.names, title, dense=True)

        # a filter string or transformation function
        if self.transf_f is not None:
            sol.apply_transf(self.transf_f, self.transf_f.names)
        return sol

def scan(model, plan,
         tf=1.0,
         npoints=500,
         t0=0.0,
         initial=None,
         times=None,
         outputs=None,
         titles=None,
         changing_pars=None):

    """Wrapper around ModelSolver."""

    plan = dict(plan)
    names = list(plan) # gets the keys
    # zip, terminating on the shortestsequence
    scan_values = list(zip(*(plan.values())))

    if titles is None:
        titles = []
        for run_values in scan_values:
            pairs = ['%s = %g'%(n, v) for (n, v) in zip(names, run_values)]
            titles.append(', '.join(pairs))

    ms = ModelSolver(model, tf=tf, npoints=npoints, t0=t0,
                     initial=initial, times=times, outputs=outputs,
                     changing_pars=names)

    s = Solutions()
    for title, run_values in zip(titles, scan_values):
        s += ms.solve(title=title, par_values=run_values)

    return s


def test():
    #import time
    from stimator import read_model
    m1_text = """
    title a simple 2 step system
    v1: A -> B, rate = V / (Km1 + A), V = 1, Km = 1
    
    v2: B ->  , rate = V * c2 * B**3
    
    V  = sqrt(4.0)
    Km1 = 1
    c2 = 0.2
    
    find c2 in [0, 1.2]
    
    init: B = 0.4, A = 1
    
    -> vin = 2 * A * v1.Km
    ~ t1 = A + B + vin
    ~ t2 = v1.V * A * step(t, 1.0)
    # ~ t3 = v1.V * A * max(t, 1.0)
    """

    m = read_model(m1_text)

    print(m1_text)

    print('\n********** Testing stoichiometry matrix ********************')
    print('Stoichiometry matrix:')
    N = genStoichiometryMatrix(m)
    print('  ', '  '.join([v.name for v in m.reactions]))
    for i, x in enumerate(m.varnames):
        print(x, N[i, :])
    print()
    print('********** Testing state2array()****************************')
    print('state2array(m):')
    v = init2array(m)
    print(v, 'of type', type(v))
    print()
    print('********** Testing rate and dXdt strings *******************')
    print('rates_strings(fully_qualified = False): ---')
    print(rates_strings(m, fully_qualified=False))
    print('\nrates_strings(): -------------------------')
    print(rates_strings(m))
    print('\ndXdt_strings(): --------------------------')
    dxdt_strs = dXdt_strings(m)
    for x in m.varnames:
        print('(d%s/dt) ='%(x), dxdt_strs[x])

    print('********** Testing differentiation of strings *******************')
    print('\n---------- Testing _gen_canonical_symbmap(m) --------------')
    symbols = _gen_canonical_symbmap(m)
    symbmap, _ = symbols['s_table'], symbols['sympy_s_table']
    print('symbmap')
    for k in symbmap:
        print('{:8} --> {}'.format(k, symbmap[k]))
    print('\n---------- Differentiation --------------')

    dxdt_strs = dXdt_strings(m)
    for x in m.varnames:
        expr = dxdt_strs[x]
        print('\nexpression =', expr)
        for x in m.varnames:
            dexpr = _differentiate_expr(expr, x, symbols)
            print('d / d', x, '=', dexpr)
        print('---')
        for p in m.parameters:
            dexpr = _differentiate_expr(expr, p.name, symbols)
            print('d / d', p.name, '=', dexpr)

    print('\nJacobian_strings(): -------------------------')
    vnames = m.varnames
    for i, vec in enumerate(Jacobian_strings(m)):
        for j, dxdx in enumerate(vec):
            print('(d d%s/dt / d %s) ='%(vnames[i], vnames[j]), dxdx)
    print('\ndfdp_strings(m, parnames): ------------------')
    parnames = "c2 v1.V".split()
    print('parnames = {}\n'.format(parnames))
    vnames = m.varnames
    for i, vec in enumerate(dfdp_strings(m, parnames)):
        for j, dxdx in enumerate(vec):
            print('(d d%s/dt / d %s) =' % (vnames[i], parnames[j]), dxdx)
    print()

    print('dfdp_strings(m, parnames): (with unknown pars)')
    parnames = "c3 v1.V".split()
    print('parnames = {}\n'.format(parnames))
    vnames = m.varnames
    for i, vec in enumerate(dfdp_strings(m, parnames)):
        for j, dxdx in enumerate(vec):
            print('(d d%s/dt / d %s) ='%(vnames[i], parnames[j]), dxdx)

    print('\n********** Testing _gen_calc_symbmap(m) *******************')
    print('_gen_calc_symbmap(m, with_uncertain = False):')
    # print(_gen_calc_symbmap(m))
    for k, v in _gen_calc_symbmap(m).items():
        print('{:8} --> {}'.format(k, v))

    print('\n_gen_calc_symbmap(m, with_uncertain = True):')
    # print(_gen_calc_symbmap(m, with_uncertain=True))
    for k, v in _gen_calc_symbmap(m, with_uncertain=True).items():
        print('{:8} --> {}'.format(k, v))
    
    print('\n********** Testing rateCalcString **************************')
    symbmap = _gen_calc_symbmap(m, with_uncertain=False)
    symbmap2 = _gen_calc_symbmap(m, with_uncertain=True)
    for v in (m.reactions.v1,
              m.reactions.v2,
              m.transformations.t1,
              m.transformations.t2,
              m.input_variables.vin):
        vstr = v(fully_qualified=True)
        print('calcstring for %s = %s\n   '% (v.name, vstr), rateCalcString(vstr, symbmap))
    print('calcstring for v2 with uncertain parameters:\n\t', rateCalcString(m.reactions.v2(fully_qualified=True), symbmap2))

    print('\n********** Testing rate and dXdt generating functions ******')
    print('Operating point --------------------------------')
    varvalues = [1.0, 0.4]
    pars = [0.4]
    t = 0.0

    print("t =", t)
    print('variables:')
    print(dict((n, value) for n, value in zip(m.varnames, varvalues)))
    print('parameters:')
    print(dict((p.name, p) for p in m.parameters))

    print('\n---- rates using all_rates_func(m) -------------------------')
    func = all_rates_func(m)
    ivs, vs, ts = func(varvalues, t)
    frmtstr = "%s = %-25s = %s"
    for v, r in zip(m.reactions, vs):
        print(frmtstr % (v.name, v(fully_qualified=True), r))
    for v, r in zip(m.transformations, ts):
        print(frmtstr % (v.name, v(fully_qualified=True), r))
    for v, r in zip(m.input_variables, ivs):
        print(frmtstr % (v.name, v(fully_qualified=True), r))

    print('---- same, at t = 2.0 --')
    ivs, vs, ts = func(varvalues, 2.0)
    for v, r in zip(m.reactions, vs):
        print(frmtstr % (v.name, v(fully_qualified=True), r))
    for v, r in zip(m.transformations, ts):
        print(frmtstr % (v.name, v(fully_qualified=True), r))
    for v, r in zip(m.input_variables, ivs):
        print(frmtstr % (v.name, v(fully_qualified=True), r))

    print('\n********** Testing add_dSdt_to_model() ***************')
    print('------ in original model')
    print('variables')
    vnames = m.varnames
    dxdtstrs = dXdt_strings(m)
    for x in vnames:
        print(x, m.get_init(x))
        print('   d {} / dt = {}'.format(x, dxdtstrs[x]))
    print('---------------------')
    pars = "Km2 v1.V init.B".split()
    print('pars =', pars)
    Snames = add_dSdt_to_model(m, pars)
    print('Snames = \n', Snames)
    print('------ in augmented model')
    print('variables')
    vnames = m.varnames
    dxdtstrs = dXdt_strings(m)
    for x in vnames:
        print(x, m.get_init(x))
        print('   d {} / dt = {}'.format(x, dxdtstrs[x]))

    print('---------------- EXAMPLE 1 ------------------')
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
    print(mtext)

    m1 = read_model(mtext)

    solution1 = solve(m1, tf=50, title='two enzymes, use !! A C ~')
    solution1a = solve(m1, tf=50, outputs='A B C sum'.split(),
                       title='explicit outputs=[A B C sum]')
    solution1v = solve(m1, tf=100, outputs='>>',
                       title='outputs=">>"')

    print('--- Last time point ----')
    print('At t =', solution1.t[-1])
    #print solution1.last
    for x in solution1.last:
        print("%-8s= %f" % (x, solution1.last[x]))

    # print 'END of EXAMPLES 1'
    # t1 = time.time()
    # print 'took', t1 - t0

    print('---------------- EXAMPLE 3 ------------------')
    m3 = read_model(models.ca.text)

    print(models.ca.text)
    ms = ModelSolver(m3, tf=8.0, npoints=2000)
    solution3 = ms.solve()
##     solution3 = solve(m3, tf = 8.0, npoints = 2000)

    # print 'END of EXAMPLES 3'
    # t3 = time.time()
    # print 'took', t3 - t1


    print('---------------- EXAMPLE 4 ------------------')
    m4 = read_model(models.rossler.text)

    print(m4)

    solution4 = solve(m4, tf=100.0, npoints=2000,
                      outputs="x1 x2 x3".split())
    solution4b = solve(m4, tf=100.0, npoints=2000, outputs="~",
                       title='Rossler, outputs="~"')

    def transformation(vars, t):
        if t > 40.0:
            return (vars[0]-5.0, vars[1], vars[2])
        else:
            return (-5.0, vars[1], vars[2])

    solution4.apply_transf(transformation,
                           new_title='Rossler, after a transformation')

    # print 'END of EXAMPLES 4'
    # t4 = time.time()
    # print 'took', t4 - t3

    # savingfile = open('examples/analysis.png', 'w+b')
    # savingfile = 'examples/analysis.png'
    sols = Solutions([solution1, solution1a, solution1v,
                      solution3,
                      solution4b, solution4])
    sols.plot() #save2file=savingfile)

    # print 'END of plotting first 4 examples'
    # tplot = time.time()
    # print 'took', tplot - t4


    print('---------------- scanning example ------------------')
    m3 = read_model(models.ca.text)
    scans = 0.0, 0.1, 0.3, 0.5, 0.8, 1.0
    # scans_k1 = 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9

    sols2 = scan(m3, {'B': scans}, tf=10.0)
    # print 'END of SCANNING EXAMPLE'
    # tscancomp = time.time()
    # print 'took', tscancomp - tplot

    sols2.plot(legend=True, ynormalize=True, group=['Ca'], fig_size=(10, 6))

    # print 'END of PLOTTING SCANNING EXAMPLE'
    # tscan = time.time()
    # print 'took', tscan - tscancomp



    print('---------------- stairway example ------------------')
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
    #print(mtext)

    mstair = read_model(mtext)

    solstairs = solve(mstair, tf=300, title='stairway')
    # print 'END of STAIRWAY EXAMPLE'
    # tstairway = time.time()
    # print 'took', tstairway - tscan

    solstairs.plot(fig_size=(9, 6), show=True)

    # print 'END of STAIRWAY PLOTTING'
    # tstairwayplot = time.time()
    # print 'took', tstairwayplot - tstairway


if __name__ == "__main__":
    test()
