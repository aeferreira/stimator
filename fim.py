#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from numpy import *
from model import *
from analysis import *

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
    tempOUT = []
    
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
            #~ tempOUT.append((Smatrix[i][j], _dres))
    

def test():
    
    #print '---------------- EXAMPLE 3 ------------------'
    m = Model("Calcium Spikes")
    m.v0 = " -> Ca", 1
    m.v1 = react(" -> Ca", rate = "B*k1")
    m.k1 = 7.3
    m.B  = 0.4
    m.export = " Ca -> ", 10
    m.leak   = "CaComp -> Ca", "kleak * CaComp"
    m.kleak  = 1.0
    m.v2     = "Ca -> CaComp", "65 * Ca**2 / (1+Ca**2)"
    m.v3     = "CaComp -> Ca", "V3*CaComp**2/(CaComp**2+K3) * Ca**4/(Ca**4 + 0.6561)"
    m.K3     = 4
    m.V3     = 500
    m.init   = state(Ca = 0.1, CaComp = 0.63655)

    print m

    #~ solution3 = solve(m3, tf = 8.0, npoints = 2000)
        
    #~ print '\nJacobian ----------------------------------------'
    #~ dxdtstrings = m.dXdt_strings()
    #~ for d in dxdtstrings:
        #~ print d
    #~ print
    #~ J = m.Jacobian_strings()
    #~ nvars = len(m.variables)
    #~ for i in range(nvars):
        #~ for j in range(nvars):
            #~ print (m.variables[i].name, m.variables[j].name),
            #~ print '\t', J[i][j]
    
    #~ print '\nJacobian evaluation -----------------------------'
    #~ J = m.getJacobian()
    #~ t = 0.0
    #~ x = array([0.1, 0.63655])
    #~ print 'with t =', t
    #~ print 'and x =', x
    #~ print 'Jacobian:'
    #~ res = J(x,t)
    #~ print res

    #~ print '\nMatrix dfdp -------------------------------------'
    #~ dxdtstrings = m.dXdt_strings()
    #~ for d in dxdtstrings:
        #~ print d
    #~ print
    #~ pars = ["B", "k1", "V3"]
    #~ print 'parameters are', pars
    #~ J = m.dfdp_strings(pars)
    #~ nvars = len(m.variables)
    #~ for i in range(nvars):
        #~ for j in range(len(pars)):
            #~ print (m.variables[i].name, pars[j]),
            #~ print '\t', J[i][j]
    
    print '\nAdding sensitivity ODEs -------------------------'
    pars = ["B", "k1", "K3"]
    add_dSdt_to_model(m, pars)
    print m
    
    plot_vars = solve(m, tf = 3.0, npoints = 900, title = 'X', outputs = "Ca CaComp")
    plot_dXdB = solve(m, tf = 3.0, npoints = 2700, title = 'dX/dB', outputs = 'd_Ca_d_B d_CaComp_d_B') #TODO: solutions should be 'clonable'.
    plot_dXdk1 = solve(m, tf = 3.0, npoints = 2700, title = 'dX/dk1', outputs = 'd_Ca_d_k1 d_CaComp_d_k1') #TODO: solutions should be 'clonable'.
    plot_dXdK3 = solve(m, tf = 3.0, npoints = 2700, title = 'dX/dK3', outputs = 'd_Ca_d_K3 d_CaComp_d_K3') #TODO: solutions should be 'clonable'.

    plot([plot_vars,plot_dXdB,plot_dXdk1,plot_dXdK3])
    


if __name__ == "__main__":
    test()

