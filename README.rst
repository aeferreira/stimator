*S-timator*: analysis of ODE models with focus on model selection and parameter estimation.
===========================================================================================

*S-timator* is a Python library to analyse ODE-based models
(also known as *dynamic* or *kinetic* models). These models are often found
in many scientific fields, particularly in Physics, Chemistry, Biology and
Engineering.

Features include:

- **A mini language used to describe models**: models can be input as plain text 
  following a very simple and human-readable language.
- **Basic analysis**: numerical solution of ODE's, parameter scanning.
- **Parameter estimation** and **model selection**: given experimental data in
  the form of time series and constrains on model operating ranges,
  built-in numerical optimizers can find parameter values and assist you in the
  experimental design for model selection.

*S-timator* is in an alpha stage: many new features will be available soon.

Requirements
------------

*S-timator* supports Python versions 2.6 and up, but support of 3.x is
coming soon.

*S-timator* depends on the "scientific python stack". The **mandatory**
requirements for *S-timator* are the following libraries:

- ``Python (2.6 or 2.7)``
- ``numpy``
- ``scipy``
- ``matplotlib``
- ``pip``

One of the following "scientific python" distributions is recommended, **as they all provide 
an easy installation of all requirements**:

- `Anaconda <https://store.continuum.io/cshop/anaconda/>`_ (or `Miniconda <http://conda.pydata.org/miniconda.html>`_ followed by the necessary ``conda install``'s)
- `Python (x,y) <https://code.google.com/p/pythonxy/>`_
- `Enthought Canopy <https://www.enthought.com/products/canopy/>`_

The installation of these Python libraries is optional, but strongly recommended:

- ``sympy``: necessary to compute dynamic sensitivities, error estimates of
  parameters and other symbolic computations.
- ``IPython`` and all its dependencies: some *S-timator* examples are provided
  as IPython notebooks.
- ``wxPython``: although *S-timator* is a python library meant to be used for scripting or in
  IPython *literate programming* interface, a simple GUI is included. This interface
  requires wxPython.


Installation
------------

After installing the required libraries, (``Python``, ``numpy``, ``scipy``,
``matplotlib`` and ``pip``) the easiest way to install *S-timator* is
with ``pip``::

    $ pip install stimator

The classical way also works, but is not recommended::
    
    $ python setup.py install

Basic use: solution of ODE models
---------------------------------

This is a warm-up example that illustrates model description, ODE numerical 
solving and plotting:

.. code:: python

    from stimator import read_model, solve

    mdl = """# Example file for S-timator
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
    !! x1 x2"""

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

    mdl = """# Example file for S-timator
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

Model selection (experimental design)
-------------------------------------

One of the examples included in *S-timator* solves an experimental design problem: 
finding a feasible set of experimental conditions that lead to the clear selection between 2 models.

Run example ``glyoxalase_discrim_2m.py``.


Summary of road map
-------------------

- Improve documentation
- I/O to other model description formats (SBML, etc)

