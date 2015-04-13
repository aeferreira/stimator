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
- ``pandas``
- ``seaborn``

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
``matplotlib``, ``pandas`` and ``pip``) the easiest way to install *S-timator* is
with ``pip``::

    $ pip install stimator

The classical way also works, but is not recommended::
    
    $ python setup.py install
