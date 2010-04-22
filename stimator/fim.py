#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2010 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
from model import *
from analysis import *
import timecourse

def add_dSdt_to_model(m, pars):
    """Add sensitivity ODEs to model, according to formula:
    
    dS/dt = df/dx * S + df/dp
    
    m is a model object
    pars are a list of parameter names
    """
    #Change par names to legal identifiers
    newpars = []
    init_of = []
    for p in pars:
        if '.' in p:
            init_of.append(p.split('.')[1])
            newpars.append(p.replace('.','_'))
        else:
            init_of.append(None)
            newpars.append(p)
    pars = newpars

    J = m.Jacobian_strings()
    dfdpstrs = m.dfdp_strings(pars)
            
    #~ print '************************'
    #~ print 'PARAMETERS'
    #~ print pars
    #~ print '************************'
    #~ print 'dfdp_strings(pars)'
    #~ for i,x in enumerate(m.variables):
        #~ for j,p in enumerate(pars):
            #~ print "d(%s, %s) = %s"%(x.name, p, dfdpstrs[i][j])
    #~ print '************************'
    nvars = len(J)
    npars = len(pars)
    
    try:
        import sympy
    except:
        print 'ERROR: sympy module must be installed to generate sensitivity strings'
        raise
    _symbs = {}
    for x in m.variables:
        _symbs[x.name] = sympy.Symbol(str(x.name))
    for p in m.parameters:
        _symbs[p.name] = sympy.Symbol(str(p.name))
    for f in m.forcing:
        _symbs[f.name] = sympy.Symbol(str(f.name))
    for p in pars:
        if not _symbs.has_key(p):
            _symbs[p] = sympy.Symbol(str(p))

    #create symbols for sensitivities
    Smatrix = []
    for i in range(nvars):
        Smatrix.append([])
        for j in range(npars):
            Sname = "d_%s_d_%s" % (m.variables[i].name, pars[j])
            _symbs[Sname] = sympy.Symbol(str(Sname))
            Smatrix[i].append(Sname)
    #~ print '************************'
    #~ print 'SYMBOLS'
    #~ for k,v in _symbs.items():
        #~ print "%-18s:"%k, v
    #~ print '************************'

    #print '************************'
    #print "SENSITIVITIES ODE's"
    for i in range(nvars):
        vname = m.variables[i].name
        for j in range(npars):
            #compute string for dS/dt
            if init_of[j] is None:
                resstr = dfdpstrs[i][j]
            else:
                resstr = ''
            for k in range(nvars):
                resstr = resstr+ "+(%s)*(%s)"%(J[i][k], Smatrix[k][j])
            _res = eval(resstr, None, _symbs)
            _dres = str(_res)
            if _dres == '0':
                _dres = '0.0'
            Smatrix[i][j] = str(Smatrix[i][j])
            #print Smatrix[i][j], '=', _dres
            setattr(m, Smatrix[i][j], variable(_dres))
            if init_of[j] is None:
                setattr(m.init, Smatrix[i][j], 0.0)
                #print Smatrix[i][j], '(0) =', 0.0
            else:
                if init_of[j] == vname:
                    setattr(m.init, Smatrix[i][j], 1.0)
                    #print Smatrix[i][j], '(0) =', 1.0
                else:
                    setattr(m.init, Smatrix[i][j], 0.0)
                    #print Smatrix[i][j], '(0) =', 0.0
    #print '************************'

def __computeNormalizedFIM(model, pars, vars, timecoursecollection, expCOV):
    
    """Computes FIM normalized by parameter values.
    
    model is a model object.
    
    pars is a dicionary of parameter names as keys and parameter values as values.
        names of the form init.<name> are possible.
        
    vars is a list of names of variables to be considered in the timecourses
    
    timecoursecollection is a Solutions or TimeCourseCollection object 
       (initial values are read from these)
    
    expCOV is an experimental variance-covariance matrix of variables.
        NOTE: only a constant error for all the variables is now possible:
        set expCOV to a float
    
    """
    m  = model.clone()
    
    #ensure m has init attr
    inits = {}
    for x in m.variables:
        inits[str(x.name)] = 0.0
    setattr(m, 'init', state(**inits))
    
    for p in pars:
        setattr(m, p, pars[p])  #TODO: verify existence of each p
    # Adding sensitivity ODEs
    npars = len(pars)
    nvars = len(vars)
    nsens = npars * nvars
    add_dSdt_to_model(m, pars.keys())

    timecoursedata = timecoursecollection
    
    #scale = float(max([ (s.t[-1]-s.t[0]) for s in timecoursedata]))
    
    sols = timecourse.Solutions()
    for tc in timecoursedata:
        #set init and solve
        for n,v in tc.state_at(0.0): #TODO: may not start at zero
            setattr(m.init, n, v)
        sols += solve(m, tf = tc.t[-1])#, scale = scale)
    
    
    #compute P and MVINV

    # P is the diagonal matrix of parameter values
    P = matrix(diag([pars[p] for p in pars]))

    # MVINV is the inverse of measurement COV matrix
    # RIGHT NOW, IT ONLY WORKS FOR A VECTOR OF CONSTANT VALUES
    # if COV is a scalar, transform in vector of constants
    if isinstance(expCOV, float) or isinstance(expCOV, int): #constant for all variables
        error_x = array([expCOV for i in range(nvars)], dtype=float)
    else:
        error_x = array(expCOV, dtype=float)
    errorvec = error_x**2
    MV = matrix(diag(errorvec))
    MVINV = linalg.inv(MV)
    
    #compute indexes of sensitivities
    # search pattern d_<var name>_d_<parname> in variable names
    indexes = []
    for vname in vars:
        for i,y in enumerate(m.variables):
            dnames = y.name.split('_')
            if len(dnames) < 3: continue
            if dnames[1] == vname and dnames[0] == 'd' and dnames[2] == 'd':
                indexes.append(i)
    indexes = array(indexes)

    tcFIM = []
    for sol in sols:
        #timestep (assumed constant)
        h = (sol.t[1] - sol.t[0]) #*scale
        #compute integral of ST * MVINV * S
        ntimes = len(sol.t)
        FIM = zeros((npars,npars))
        for i in range(ntimes):
            #compute S matrix
            svec = sol.data[indexes , i]
            S = matrix(svec.reshape((nvars, npars)))
            S = S * P # scale with par values
            ST = matrix(S.T)

            #compute contribution at point i and add to FIM
            FIMpoint = h * ST * MVINV * S
            FIM += FIMpoint
        tcFIM.append(FIM)
    
    sumFIM = zeros((npars,npars))  #TODO use numpy tensor and 1 dimension sum?
    for f in tcFIM:
        sumFIM += f
    
    return sumFIM, P



def __computeNormalizedFIM4SingleTC(model, pars, vars, tf, expCOV):
    
    sols = Solutions()
    
    sols += solve(model, tf)
    FIM, P = __computeNormalizedFIM(model, pars, vars, sols, expCOV)
    return FIM, P
        
def computeFIM4SingleTC(model, pars, vars, tf, COV):
    FIM, P = __computeNormalizedFIM4SingleTC(model, pars, vars, tf, COV)
    # compute inverse of P to "descale"
    PINV = linalg.inv(P)
    # compute FIM and its inverse (lower bounds for parameter COV matrix)
    realFIM = PINV * FIM * PINV
    INVFIM = linalg.inv(FIM)
    INVFIM = P * INVFIM * P
    check = dot(INVFIM, realFIM)
    #TODO: MAKE THIS CHECK USEFUL
    
    return realFIM, INVFIM

def computeFIM(model, pars, vars, TCs, COV):
    FIM, P = __computeNormalizedFIM(model, pars, vars, TCs, COV)
    # compute inverse of P to "descale"
    PINV = linalg.inv(P)
    # compute FIM and its inverse (lower bounds for parameter COV matrix)
    realFIM = PINV * FIM * PINV
    INVFIM = linalg.inv(FIM)
    INVFIM = P * INVFIM * P
    check = dot(INVFIM, realFIM)
    #TODO: MAKE THIS CHECK USEFUL
    return realFIM, INVFIM

def test():
    
    m = Model("Glyoxalase system")
    m.glo1 = react("HTA -> SDLTSH", rate = "V1*HTA/(Km1 + HTA)")
    m.glo2 = react("SDLTSH -> "   , rate = "V2*SDLTSH/(Km2 + SDLTSH)")
    m.V1  = 2.57594e-05
    m.Km1 = 0.252531
    m.V2  = 2.23416e-05
    m.Km2 = 0.0980973
    m.init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)
    
    pars = "V1 Km1".split()
    parvalues = [getattr(m, p) for p in pars]
    
    sols = timecourse.Solutions()
    sols += solve(m, tf = 4030.0) 
    
    pars = "V1 Km1".split()
    parvalues = [getattr(m, p) for p in pars]
    parsdict = dict (zip(pars, parvalues))
    
    errors = (0.01,0.001)

    print '\n------------------------------------------------'
    print 'Glyoxalase model, 1 timecourse, parameters V1 and Km1'
    print 'Timecourse with HTA and SDLTSH'
    FIM1, invFIM1 = computeFIM(m, parsdict, "HTA SDLTSH".split(), sols, errors)
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    print '\n------------------------------------------------'
    print 'Glyoxalase model, 1 timecourse, parameters V1 and Km1'
    print 'Timecourse with SDLTSH only'
    FIM1, invFIM1 = computeFIM(m, parsdict, "SDLTSH".split(), sols, errors[1])
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    

    print '\n------------------------------------------------'
    print 'Glyoxalase model, 2 timecourses, parameters V1 and Km1'
    print 'Timecourses with HTA and SDLTSH'

    #generate 2nd timecourse
    m.init.SDLTSH = 0.001246154
    m.init.HTA = 0.2688
    sols += solve(m, tf = 5190.0)
    
    FIM1, invFIM1 = computeFIM(m, parsdict, "HTA SDLTSH".split(), sols, errors)
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    
    print '\n------------------------------------------------'
    print 'Glyoxalase model, 2 timecourses, parameters V1 and Km1'
    print 'Timecourses with SDLTSH only'
    
    FIM1, invFIM1 = computeFIM(m, parsdict, ["SDLTSH"], sols, errors[1])
    print '\nParameters ---------------------------'
    for i,p in enumerate(parsdict.keys()):
        print "%7s = %.5e +- %.5e" %(p, parsdict[p], invFIM1[i,i]**0.5)    
    

if __name__ == "__main__":
    test()

