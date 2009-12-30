#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
from model import *
from analysis import *
import timecourse

def add_dSdt_to_model(m, pars):
    """Add sensitivity ODS to model, according to formula:
    
    dS/dt = df/dx * S + df/dp
    
    m is a model object
    pars are a list of parameter names
    """
    J = m.Jacobian_strings()
    dfdpstrs = m.dfdp_strings(pars)
    nvars = len(J)
    npars = len(pars)
    
    try:
        import sympy
    except:
        print 'ERROR: sympy module must be installed to generate sensitivity strings'
        raise
    _symbs = {}
    for x in m.variables:
        _symbs[x.name] = sympy.Symbol(x.name)
    for p in m.parameters:
        _symbs[p.name] = sympy.Symbol(p.name)
    for f in m.forcing:
        _symbs[f.name] = sympy.Symbol(f.name)
    
    #create symbols for sensitivities
    Smatrix = []
    for i in range(nvars):
        Smatrix.append([])
        for j in range(npars):
            Sname = "d_%s_d_%s" % (m.variables[i].name, pars[j])
            _symbs[Sname] = sympy.Symbol(Sname)
            Smatrix[i].append(Sname)

    for i in range(nvars):
        for j in range(npars):
            #compute string for dx/dt
            resstr = dfdpstrs[i][j]
            for k in range(nvars):
                resstr = resstr+ "+(%s)*(%s)"%(J[i][k], Smatrix[k][j])
            #~ print Smatrix[i][j]
            #~ print resstr
            #~ print type(resstr)
            _res = eval(resstr, None, _symbs)
            _dres = str(_res)
            if _dres == '0':
                _dres = '0.0'
            setattr(m, Smatrix[i][j], variable(_dres))
            setattr(m.init, Smatrix[i][j], 0.0)

def __computeScaledFIM(model, pars, vars, timecoursecollection, expCOV):
    
    """Computes FIM normalized by parameter values.
    
    model is a model object.
    
    pars is a dicionary of parameter names as keys and parameter values as values.
        names of the form init.<name> are possible.
        
    vars is a list of names of variables to be considered in the timecourses
    
    timecoursecollection is a TimeCourseCollection object 
       (initial values and parameters are read from these)
    
    expCOV is an experimental variance-covariance matrix of variables.
        NOTE: only a constant error for all the variables is now possible:
        set expCOV to a float
    
    """
    m  = model.clone()
    
    for p in pars:
        setattr(m, p, pars[p])  #TODO: verify existence of each p

    # Adding sensitivity ODEs
    npars = len(pars)
    nvars = len(vars)
    nsens = npars * nvars
    add_dSdt_to_model(m, pars.keys())
    
    #sol = SolutionTimeCourse (d[:,0].T, d[:,1:].T, h[1:])
    
    if (isinstance(timecoursecollection, timecourse.TimeCourseCollection)):
        # change to SolutionTimeCourse interface
        timecoursedata = Solutions()
        for i,d in enumerate(timecoursecollection.data):
            timecoursedata += SolutionTimeCourse (d[:,0].T, d[:,1:].T, timecoursecollection.headers[i][1:])
            
    else:
        timecoursedata = timecoursecollection

    scale = float(max([ (s.t[-1]-s.t[0]) for s in timecoursedata]))
    
    calcDerivs = model.getdXdt(scale=scale, with_uncertain=False)
    salg=integrate._odepack.odeint
    names = [x.name for x in m.variables]

    X0 = []
    times = []
    varindexes = []
    for d in timecoursedata:
        y0 = copy(d.data[:, 0]) # 1st col are initial values
        X0.append(y0)
        
        t  = d.t 
        t0 = t[0]
        tts = (t-t0)/scale+t0  # scale time points
        times.append(tts)

    #Computes solution for timecourse i, given parameters trial.
    sols = Solutions()
    for i in range(len(timecoursedata)):
        y0 = copy(X0[i])
        # fill uncertain initial values
            
        t  = copy(times[i])

        #~ Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
        output = salg(calcDerivs, y0, t, (), None, 0, -1, -1, 0, None, 
                        None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
        #if output[-1] < 0: return None
        Y = output[0]
        title = "timecourse %d"%i
        s = SolutionTimeCourse (t, Y.T, names, title)
        sols+=s
    
    plot( sols )
    return
    
    #compute P
    # P is the diagonal matrix of parameter values
    parvalues = [getattr(m, p) for p in pars]
    P = matrix(diag(parvalues))
    
    #compute MVINV
    # MVINV is the inverse of measurement COV matrix
    # RIGHT NOW, IT ONLY WORKS FOR A SINGLE CONSTANT VALUE
    error_x = expCOV
    errorvar = error_x**2
    errorvec = [errorvar for i in range(nvars)]
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
        h = (sol.t[1] - sol.t[0]) *scale
        #compute integral of ST * MVINV * S
        ntimes = sol.data.shape[1]
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
    
    sumFIM = zeros((npars,npars))
    for f in tcFIM:
        sumFIM += f
    
    return sumFIM, P



def __computeNormalizedFIM(model, pars, vars, tf, expCOV):
    m = model.clone()

    # Adding sensitivity ODEs
    npars = len(pars)
    nvars = len(vars)
    nsens = npars * nvars
    add_dSdt_to_model(m, pars)
    
    #solve complete model
    sol = solve(m, tf = tf)

    #timestep (assumed constant)
    h = (sol.t[1] - sol.t[0])
    
    #compute P
    # P is the diagonal matrix of parameter values
    parvalues = [getattr(m, p) for p in pars]
    P = matrix(diag(parvalues))
    
    #compute MVINV
    # MVINV is the inverse of measurement COV matrix
    # RIGHT NOW, IT ONLY WORKS FOR A SINGLE CONSTANT VALUE
    error_x = expCOV
    errorvar = error_x**2
    errorvec = [errorvar for i in range(nvars)]
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
    
    #compute integral of ST * MVINV * S
    ntimes = sol.data.shape[1]
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
    
    return FIM, P
    
def computeFIM(model, pars, vars, tf, expCOV):
    FIM, P = __computeNormalizedFIM(model, pars, vars, tf, expCOV)
    # compute inverse of P to "descale"
    PINV = linalg.inv(P)
    # compute FIM and its inverse (lower bounds for parmeter COV matrix)
    realFIM = PINV * FIM * PINV
    INVFIM = linalg.inv(FIM)
    INVFIM = P * INVFIM * P
    check = dot(INVFIM, realFIM)
    #TODO: MAKE THIS CHECK USEFUL
    
    return realFIM, INVFIM

def computeFIM4TCs(model, pars, vars, timecoursecollection, expCOV):
    FIM, P = __computeScaledFIM(model, pars, vars, timecoursecollection, expCOV)
    # compute inverse of P to "descale"
    PINV = linalg.inv(P)
    # compute FIM and its inverse (lower bounds for parmeter COV matrix)
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


    print '\n--------- Example 1 ------------------'
    print 'Glyoxalase model, 1 timecourse, parameters V1 and Km1'
    print 'Timecourse with HTA and SDLTSH'
    print
    FIM1, invFIM1 = computeFIM(m, pars, "HTA SDLTSH".split(), 4030.0, 0.001)
    print "FIM ="
    print FIM1
    print 'FIM-1 ='
    print invFIM1
    print '\nParameters ---------------------------'
    for i in range(len((pars))):
        print "%7s = %.5e +- %.5e" %(pars[i], parvalues[i], invFIM1[i,i]**0.5)
    print '\n--------- Example 2 ------------------'
    print 'Glyoxalase model, 1 timecourse, parameters V1 and Km1'
    print 'Timecourse with SDLTSH only'
    print
    FIM1, invFIM1 = computeFIM(m, pars, "SDLTSH".split(), 4030.0, 0.001)
    print "FIM ="
    print FIM1
    print 'FIM-1 ='
    print invFIM1
    print '\nParameters ---------------------------'
    for i in range(len((pars))):
        print "%7s = %.5e +- %.5e" %(pars[i], parvalues[i], invFIM1[i,i]**0.5)
    
    print '------------------------------------------------'
    print '|    Using new interface                       |'
    print '------------------------------------------------'
    print '\n--------- Example 1 ------------------'
    print 'Glyoxalase model, 1 timecourse, parameters V1 and Km1'
    print 'Timecourse with HTA and SDLTSH'
    print
    
    sols = Solutions()
    sols += solve(m, tf = 4030.0) 

    pars = "V1 Km1".split()
    parvalues = [getattr(m, p) for p in pars]
    pars = dict (zip(pars, parvalues))

    FIM1, invFIM1 = computeFIM4TCs(m, pars, "HTA SDLTSH".split(), sols, 0.001)
    print "FIM ="
    print FIM1
    print 'FIM-1 ='
    print invFIM1
    print '\nParameters ---------------------------'
    for i in range(len((pars))):
        print "%7s = %.5e +- %.5e" %(pars[i], parvalues[i], invFIM1[i,i]**0.5)    
    


if __name__ == "__main__":
    test()

