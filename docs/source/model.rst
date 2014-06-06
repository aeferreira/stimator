Primer for documentation on ``model.py``
****************************************

``model.py`` includes class ``Model`` and utility functions for describing a kinetic network with possible uncertain parameters.

1- Building a model
===================

Example::

    from model import *

    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")

    m.t1 = transf("A*4 + C")
    m.t2 = transf("sqrt(2*A)")

    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.V3.uncertainty(0.1, 1.0) #the same as last line
    m.Km3 = 4

    m.init       = state(A = 1.0, C = 1)
    m.afterwards = state(A = 1.0, C = 2)
    m.afterwards.C.uncertainty(1,3)

``Model()`` constructor string argument is a title for the model.

All model data is inserted through assignment (using the mighty dot). 

react() creates reactions. First argument is a string with stoichiometry. 
rate argument (second argument) is a string with rate law. 
If an int or a float, mass-action kinetics is assumed, 
and the value is the rate constant.

Model variables are guessed from the stoichiometry of reactions.

``transf()`` creates time varying expressions. Single argument (a string) is the 
expression to be computed. It can depend on variables and parameters.

Parameters are created by assignment to constant ``float``'s or ``int``'s.

If a variable (guessed by stoichiometry) is assigned to a constant it becomes 
a parameter (an 'external' or 'independent' variable). In the example, 
B is such a variable.

Parameters to be estimated are flagged with bounds for uncertainty by 
assigning them to a pair of floats (list or tuple) or by using function 
uncertainty(min, max). In the example, V3, was flagged as "uncertain" within 
0.1 and 1.0.

state function creates a state. A state is a given set of values for the 
variables of the model. Arguments of state are name = value. Only names of 
variables can appear. Missing names default to 0.0.
The main use for states is the indication of initial values in simulations, 
but they are really a multidimensional point for the dynamic system underlying 
the model with several potential uses.

In the example, there are two states, 'init' and 'afterwards'. 'init' should 
be used as the name for the default initial value state.

In the example, the variable C of state afterwards has been flagged to be 
estimated, by using function uncertainty.

2- printing and cloning
=======================

Models can be printed with print and cloned with function ``clone()``::

    print m

    Prints:
    My first model

    Variables: A C
    External variables: B
    v1:
      reagents: [('A', 1.0), ('B', 1.0)]
      products: [('C', 1.0)]
      rate    = 3*A*B
    v2:
      reagents: []
      products: [('A', 1.0)]
      rate    = 1.0
    v3:
      reagents: [('C', 1.0)]
      products: []
      rate    = V3 * C / (Km3 + C)
    t1:
      rate = A*4 + C
    t2:
      rate = sqrt(2*A)
    init: (A = 1.0, C = 1.0)
    afterwards: (A = 1.0, C = 2.0)
    B = 2.2
    myconstant = 4.0
    V3 = 0.5
    Km3 = 4.0
    V3 = ? (0.1, 1.0)
    afterwards.C = ? (1.0, 3.0)

    m2 = m.clone()
    print m2

Prints the same data!


3-Component retrieval
=====================

Most of the data in a Model can be retrieved by using the dot 
access (attribute referencing).

Most components have a name attribute. (Model has a title). 

Parameters and values in states behave has floats. 

States can be iterated with for. (name, value) tuples are returned.::

    print '********** Testing component retrieval *********************'
    print 'm.K3 :',m.Km3
    print 'm.K3.name :',m.Km3.name, '(a float with a name attr)'
    print 'm.init:',m.init
    print 'm.init.A :',m.init.A
    print 'iterating m.init'
    for name, x in m.init:
        print '	', name, '=', x
    print

    Prints:
    m.K3 : 4.0
    m.K3.name : Km3 (a float with a name attr)
    m.init: (A = 1.0, C = 1.0)
    m.init.A : 1.0
    iterating m.init
        A = 1.0
        C = 1.0
    print '********** Testing component reassignment *****************'
    print 'm.myconstant :',m.myconstant
    print len(m.parameters), 'parameters total'
    print 'making m.myconstant = 5.0'
    m.myconstant = 5.0
    print 'm.myconstant :',m.myconstant
    print len(m.parameters), 'parameters total'

    print 'making m.myconstant = react("A+B -> C"  , 3)'
    try:
        m.myconstant = react("A+B -> C"  , 3)
    except BadTypeComponent:
        print 'Failed! BadTypeComponent was caught.'
    print 'm.myconstant :',m.myconstant, '(still!)'
    print len(m.parameters), 'parameters total'
    print
    print 'm.V3 :', m.V3
    print 'm.V3.bounds:' , m.V3.bounds
    print 'iterating m.uncertain'
    for x in m.uncertain:
        print '	', x.name, 'in (', x.min, ',', x.max, ')'
    print len(m.uncertain), 'uncertain parameters total'
    print 'making m.V3 = [0.1, 0.2]'
    m.V3 = [0.1, 0.2]
    print 'm.V3 :', m.V3
    print 'm.V3.bounds:' ,m.V3.bounds
    print len(m.uncertain), 'uncertain parameters total'
    print 'making m.V4 = [0.1, 0.6]'
    m.V4 = [0.1, 0.6]
    print 'm.V4 :', m.V4
    print 'm.V4.bounds:' ,m.V4.bounds
    print len(m.uncertain), 'uncertain parameters total'
    print 'iterating m.uncertain'
    for x in m.uncertain:
        print '	', x.name, 'in (', x.min, ',', x.max, ')'
    print 'making m.init.A = 5.0'
    m.init.A = 5.0
    print 'iterating m.init'
    for name, x in m.init:
        print '	', name, '=', x.pprint()
    print 'flagging init.A as uncertain with   m.init.A = (0.5, 2.5)'
    m.init.A = (0.5, 2.5)
    print 'iterating m.init'
    for name, x in m.init:
        print '	', name, '=', x.pprint()
    print 'calling    m.init.A.uncertainy(0.5,3.0)'
    m.init.A.uncertainty(0.5,3.0)
    print 'iterating m.init'
    for name, x in m.init:
        print '	', name, '=', x.pprint()
    print 'calling    m.init.A.uncertainy(None)'
    m.init.A.uncertainty(None)
    print 'iterating m.init'
    for name, x in m.init:
        print '	', name, '=', x.pprint()
    print 'making m.init.A back to 1.0'
    m.init.A = 1.0
    print 'iterating m.init'
    for name, x in m.init:
        print '	', name, '=', x.pprint()
    print 
    print '********** Testing stoichiometry matrix ********************'
    print 'Stoichiometry matrix:'
    N = m.genStoichiometryMatrix()
    print '  ', '  '.join([v.name for v in m.reactions])
    for i,x in enumerate(m.variables):
        print x.name, N[i, :]
    print
    print '********** Testing rateCalcString **************************'
    print 'calcstring for v3:
    ', m.rateCalcString(m.v3.rate)
    print
    print 'calcstring for v3 with uncertain parameters:
    ', m.rateCalcString(m.v3.rate, True)
    print
    print '********** Testing rate and dXdt generating functions ******'
    print 'Operating point:'
    varvalues = [1.0, 1.0]
    pars      = [1.0]
    print 'variables  =', dict((v.name, value) for v,value in zip(m.variables, varvalues))
    print 'parameters =', dict((p.name, p)     for p in m.parameters)
    print '---- rates using Model.rates_func() -------------------------'
    vratesfunc = m.rates_func()
    vrates = vratesfunc(varvalues,0)
    for v,r in zip(m.reactions, vrates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)
    print '---- transformations using Model.transf_func() --------------'
    tratesfunc = m.transf_func()
    trates = tratesfunc(varvalues,0)
    for v,r in zip(m.transf, trates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)
    print '---- dXdt using Model.dXdt() --------------------------------'
    #f = m.getdXdt()
    dXdt = m.dXdt(varvalues,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)
    print '---- dXdt using Model.dXdt() setting uncertain parameters ---'
    print 'f = m.getdXdt(with_uncertain = True)'
    f = m.getdXdt(with_uncertain = True)
    print 'setting uncertain as', dict((v.name, value) for v,value in zip(m.uncertain, pars))
    print 'm.set_uncertain(pars)'
    m.set_uncertain(pars)
    dXdt = f(varvalues,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)
    print '---- dXdt using Model.dXdt_with(pars) ------------------------'
    print 'f = m.dXdt_with(pars)'
    f = m.dXdt_with(pars)
    dXdt   = f(varvalues,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)
    print '---- dXdt using Model.dXdt() with a state argument (m.init) --'
    print 'm.init:', m.init
    print 'making m.V3 = 1.0'
    m.V3 = 1.0
    print 'm.V3 :', m.V3
    print
    print 'f = m.dXdt'
    f = m.dXdt
    print 'dXdt = f(m.vectorize("init"),0)'
    dXdt = f(m.vectorize("init"),0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)
    print '---- same, changing state argument ---------------------------'
    m.init.A = 2.0
    print 'after m.init.A = 2.0'
    print 'm.init:', m.init
    print
    print 'f = m.dXdt'
    f = m.dXdt
    print 'dXdt = f(m.vectorize("init"),0)'
    dXdt = f(m.vectorize("init"),0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

All docs, for now. IPython based docs are comming soon.
