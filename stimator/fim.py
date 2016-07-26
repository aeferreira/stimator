"""functions to compute Fisher Information Matrix"""
from __future__ import print_function, absolute_import
from collections import OrderedDict
from numpy import array, diag, matrix, zeros, linalg, dstack, sum
from stimator.dynamics import add_dSdt_to_model, dXdt_strings, solve
from stimator.timecourse import Solutions, constError_func

SYMPY_INSTALLED = True
try:
    import sympy
except ImportError:
    print ('ERROR: sympy module must be installed to generate sensitivity strings')
    SYMPY_INSTALLED = False

def __computeNormalizedFIM(model, pars, timecoursedata, expCOV, vars=None):

    """Computes FIM normalized by parameter values.

    model is a model object.

    pars is a dicionary of parameter names as keys and parameter values as values,
    or list of (name,value) pairs
        names of the form init.<name> are possible.

    timecoursecollection is a Solutions object
       (initial values are read from these)

    expCOV is a function wich computes the
        experimental variance-covariance matrix of variables.
        function has signature f(x), where x is the numpy array of variables.
        NOTE: module expcov defines some 'stock' error functions.

    vars is a list of names of variables to be considered in the timecourses

    """
    m = model.copy()

    #ensure m has init attr
    inits = {}
    for x in m.varnames:
        inits[str(x)] = 0.0
    m.set_init(inits)
    
    convert_pars = OrderedDict()
    if isinstance(pars, dict):
        pars = pars.items()
    for p, v in pars:
        convert_pars[p] = v

    pars = convert_pars

    for n, v in pars.items():
        if '.' in n:
            alist = n.split('.')
            vname, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vname == 'init':
                m.set_init([(name, v)])
                continue
        m.setp(n, v)
    # Adding sensitivity ODEs
    npars = len(pars)
    nvars = len(vars)
    Snames = add_dSdt_to_model(m, list(pars.keys()))

    sols = Solutions()
    for tc in timecoursedata:
        #set init and solve
        m.set_init(tc.state_at(0.0))
        sols += solve(m, tf=tc.t[-1], ignore_replist=True)

    #compute P

    # P is the diagonal matrix of parameter values
    P = matrix(diag(array(list(pars.values()), dtype=float)))

    #keep indexes of variables considered
    if vars is not None:
        vnames = m.varnames
        xindexes = []
        for vname in vars:
            for i, y in enumerate(vnames):
                if y == vname:
                    xindexes.append(i)
        xindexes = array(xindexes)
        #compute indexes of sensitivities
        # search pattern d_<var name>_d_<parname> in variable names
        indexes = []
        for vname in vars:
            for Sname in Snames:
                if Sname[0] == vname:
                    indexes.append(vnames.index(Sname[2]))
        indexes = array(indexes)
        # insert indexes in sols as attributes
        for sol in sols:
            sol.xindexes = xindexes
            sol.indexes = indexes

    tcFIM = []
    for sol in sols:
        #compute integral of ST * MVINV * S
        ntimes = len(sol.t)
        FIM = zeros((npars, npars))
        for i in range(1, ntimes):
            #time step
            h = (sol.t[i]-sol.t[i-1])

            #S matrix
            svec = sol.data[sol.indexes, i]
            S = matrix(svec.reshape((nvars, npars)))
            SP = S * P # scale with par values
            SPT = matrix(SP.T)

            # MVINV is the inverse of measurement COV matrix
            xvec = sol.data[sol.xindexes, i]
            error_x = expCOV(xvec)
            MV = matrix(error_x**2)
            MVINV = linalg.inv(MV)

            #contribution at point i (rectangle's rule)
            FIMpoint = h * SPT * MVINV * SP
            FIM += FIMpoint
        tcFIM.append(FIM)

    sumFIM = sum(dstack(tcFIM), axis=2)

    return sumFIM, P

def computeFIM(model, pars, TCs, COV, vars=None):
    FIM, P = __computeNormalizedFIM(model, pars, TCs, COV, vars)
    # compute inverse of P to "descale"
    PINV = linalg.inv(P)
    # compute FIM and its inverse (lower bounds for parameter COV matrix)
    realFIM = PINV * FIM * PINV
    INVFIM = linalg.inv(FIM)
    INVFIM = P * INVFIM * P
    # check = dot(INVFIM, realFIM)
    # TODO: MAKE THIS CHECK USEFUL
    return realFIM, INVFIM

def test():
    from stimator.model import Model
    m = Model("Glyoxalase system")
    m.set_reaction('glo1', "HTA -> SDLTSH", "V*HTA/(Km1 + HTA)", pars=dict(V=2.57594e-05))
    m.set_reaction('glo2', "SDLTSH -> ", "V2*SDLTSH/(Km2 + SDLTSH)")
    m.setp('V', 2.57594e-05)
    m.setp('Km1', 0.252531)
    m.setp('V2', 2.23416e-05)
    m.setp('Km2', 0.0980973)
    m.set_init(SDLTSH=7.69231E-05, HTA=0.1357)

    pars = "glo1.V Km1".split()
    parvalues = [m.parameters.glo1.V, m.parameters.Km1]

    sols = Solutions()
    sols += solve(m, tf=4030.0)

    parsdict = dict(zip(pars, parvalues))

    errors = (0.01, 0.001)
    errors = constError_func(errors)
    errors_sdl_only = 0.001
    errors_sdl_only = constError_func(errors_sdl_only)

    divider = '\n------------------------------------------------'

    print(divider)
    glo_1tc_msg = 'Glyoxalase model, 1 timecourse, parameters {} and {}'.format
    glo_2tc_msg = 'Glyoxalase model, 2 timecourses, parameters {} and {}'.format
    print (glo_1tc_msg(pars[0], pars[1]))
    print ('Timecourse with HTA and SDLTSH')
    _, invFIM1 = computeFIM(m, parsdict, sols, errors, "HTA SDLTSH".split())

    print ('\nParameters ---------------------------')
    for i, p in enumerate(parsdict.keys()):
        print ("%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i, i]**0.5))
    print(divider)
    print (glo_1tc_msg(pars[0], pars[1]))
    print ('Timecourse with SDLTSH only')
    _, invFIM1 = computeFIM(m, parsdict, sols, errors_sdl_only, "SDLTSH".split())
    print ('\nParameters ---------------------------')
    for i, p in enumerate(parsdict.keys()):
        print ("%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i, i]**0.5))


    print(divider)
    print(glo_2tc_msg(pars[0], pars[1]))
    print('Timecourses with HTA and SDLTSH')

    #generate 2nd timecourse
    m.init.SDLTSH = 0.001246154
    m.init.HTA = 0.2688
    sols += solve(m, tf=5190.0)

    _, invFIM1 = computeFIM(m, parsdict, sols, errors, "HTA SDLTSH".split())
    print ('\nParameters ---------------------------')
    for i, p in enumerate(parsdict.keys()):
        print ("%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i, i]**0.5))

    print(divider)
    print(glo_2tc_msg(pars[0], pars[1]))
    print ('Timecourses with SDLTSH only')

    _, invFIM1 = computeFIM(m, parsdict, sols, errors_sdl_only, ["SDLTSH"])
    print ('\nParameters ---------------------------')
    for i, p in enumerate(parsdict.keys()):
        print ("%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i, i]**0.5))

if __name__ == "__main__":
    test()

