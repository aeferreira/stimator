.. _installing:

Python and "Scientific" Python
------------------------------

*S-timator* is a (pure) Python package that is installable
from the `Python Package Index <https://pypi.python.org/pypi>`_.

However, Python relies heavily on the "Scientific Python ecossystem", a set of
Python libraries that brings high-performance scientific computing to the Python
programming language.

In practical terms, this means that the installation of a standard Python
distribution will not be enough to use S-timator. Instead, either the
mandatory dependencies are installed one by one, or, more conveniently, Python
is installed through a "scientific distribution".

One of the following "scientific python" distributions is recommended,
**as they all provide an easy installation of all requirements**:

- `Anaconda <https://store.continuum.io/cshop/anaconda/>`_(or `Miniconda <http://conda.pydata.org/miniconda.html>`_ followed by the necessary ``conda install``'s)
- `Python (x,y) <https://code.google.com/p/pythonxy/>`_
- `Enthought Canopy <https://www.enthought.com/products/canopy/>`_


Requirements
~~~~~~~~~~~~

*S-timator* supports Python versions 2.6 and up, but support of 3.x is
coming soon.

*S-timator* depends on a "scientific python stack". The **mandatory**
requirements for *S-timator* are the following libraries:

- ``Python (2.6 or 2.7)``
- ``numpy``
- ``scipy``
- ``matplotlib``
- ``pip``


The installation of these Python libraries is optional, but strongly recommended:

- ``sympy``: necessary to compute dynamic sensitivities, error estimates of
  parameters and other symbolic computations.
- ``IPython`` and all its dependencies: some *S-timator* examples are provided
  as IPython notebooks.
- ``wxPython``: although *S-timator* is a python library meant to be used for scripting or in
  IPython notebook *literate computing* interface, a simple GUI is included.
  This interface requires wxPython.

The `Anaconda Python Distribution, from "Continuum Analytics" <https://store.continuum.io/cshop/anaconda/>`_
is, arguably, the most convenient distribution. The full installation will provide
all S-timator requirements, except wxPython, which has to be installed after
installing Python.

From the same company, the `Miniconda <http://conda.pydata.org/miniconda.html>`_ "slim" distribution
is also an alternative, for those that worry about disk sapce. In this case,
the necessary ``conda install``'s must be run for the dependencies, after installing Python.

Installation of S-timator
~~~~~~~~~~~~~~~~~~~~~~~~~

*S-timator* is on the Python Package Index (pypi), so , after installing the
required libraries, (``Python``, ``numpy``, ``scipy``,
``matplotlib`` and ``pip``) the easiest way to install *S-timator* is
with ``pip``::

    $ pip install stimator



