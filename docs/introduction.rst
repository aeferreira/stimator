.. _introduction:

Basic use
=========


Solution of ODE models
----------------------

This is a warm-up example that illustrates model description, ODE numerical 
solving and plotting:

.. code:: python

    from stimator import read_model, solve

    mdl = """
    # Example file for S-timator
    title Example 1

    #reactions (with stoichiometry and rate)
    vin  : -> x1     , rate = k1
    v2   : x1 ->  x2 , rate = k2 * x1
    vout : x2 ->     , rate = k3 * x2

    #parameters and initial state
    k1 = 1
    k2 = 2
    k3 = 1
    init: (x1=0, x2=0)

    #filter what you want to plot
    !! x1 x2
    """

    m = read_model(mdl)

    print '========= model ========================================'
    print mdl
    print '--------------------------------------------------------'

    solve(m, tf=5.0).plot(show=True)

Parameter estimation
--------------------

Model parameter estimation, based on experimental time-course data 
(run example ``par_estimation_ex2.py``):

.. code:: python

    from stimator import read_model, readTCs, solve
    from stimator.deode import DeODEOptimizer

    mdl = """
    # Example file for S-timator
    title Example 2

    vin  : -> x1     , rate = k1
    v2   : x1 ->  x2 , rate = k2 * x1
    vout : x2 ->     , rate = k3 * x2

    init : x1=0, x2=0
    !! x2
    find k1 in [0, 2]
    find k2 in [0, 2]
    find k3 in [0, 2]

    timecourse ex2data.txt
    generations = 200   # maximum generations for GA
    genomesize = 60     # population size in GA
    """
    m1 = read_model(mdl)
    print mdl

    optSettings={'genomesize':60, 'generations':200}
    timecourses = readTCs(['ex2data.txt'], verbose=True)

    optimizer = DeODEOptimizer(m1,optSettings, timecourses)
    optimizer.run()
    
    best = optimizer.optimum
    print best.info()
    best.plot()

This produces the following output::

    -------------------------------------------------------
    file .../examples/ex2data.txt:
    11 time points, 2 variables    

    Solving Example 2...
    0   : 3.837737
    1   : 3.466418
    2   : 3.466418
    ...  (snip)
    39  : 0.426056
    refining last solution ...

    DONE!
    Too many generations with no improvement in 40 generations.
    best energy = 0.300713
    best solution: [ 0.29399228  0.47824875  0.99081065]
    Optimization took 8.948 s (00m 08.948s)

    --- PARAMETERS           -----------------------------
    k3	    0.293992 +- 0.0155329
    k2	    0.478249 +- 0.0202763
    k1	    0.990811 +- 0.0384208

    --- OPTIMIZATION         -----------------------------
    Final Score	0.300713
    generations	40
    max generations	200
    population size	60
    Exit by	Too many generations with no improvement


    --- TIME COURSES         -----------------------------
    Name		Points		Score
    ex2data.txt	11	0.300713

