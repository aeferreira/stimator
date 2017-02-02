.. _installing:

Installation: Python and "Scientific" Python
--------------------------------------------

*S-timator* is a (pure) Python package that is installable
from the `Python Package Index <https://pypi.python.org/pypi>`_.

However, *S-timator* relies heavily on the "Scientific Python ecossystem", a set of
Python libraries that brings high-performance scientific computing to the Python
programming language.

In practical terms, this means that the installation of a standard Python
distribution will not be enough to use S-timator. Instead, either the
mandatory dependencies are installed one by one, or, more conveniently, Python
is installed through a "scientific distribution".

One of the following "scientific python" distributions is recommended,
**as they all provide an easy installation of all requirements**:

- `Anaconda <https://store.continuum.io/cshop/anaconda/>`_ (or `Miniconda <http://conda.pydata.org/miniconda.html>`_ followed by the necessary ``conda install``'s)
- `Python (x,y) <https://code.google.com/p/pythonxy/>`_
- `Enthought Canopy <https://www.enthought.com/products/canopy/>`_


Requirements
~~~~~~~~~~~~

*S-timator* supports Python versions 2.7 and 3.3 or above.

*S-timator* depends on a "scientific python stack". The **mandatory**
requirements for *S-timator* are the following libraries:

- ``Python``
- ``numpy``
- ``scipy``
- ``matplotlib``
- ``pip``
- ``pandas``
- ``seaborn``


Other Python libraries which are optional, but strongly recommended:

- ``sympy``: necessary to compute dynamic sensitivities, error estimates of
  parameters and other symbolic computations.
- The ``Jupyter`` package: some *S-timator* examples are provided
  as Jupyter notebooks.

The `Anaconda Python Distribution, from "Continuum Analytics" <https://store.continuum.io/cshop/anaconda/>`_
is a convenient Python distribution that includes all of the above requirements.

If the full `Anaconda Distribution` is too heavy on disk space, from the same company, the `Miniconda <http://conda.pydata.org/miniconda.html>`_ "slim" distribution
is also an alternative. In this case, run the necessary ``conda install``'s of the requirements, after installing Miniconda:

    $ conda install numpy scipy sympy pandas seaborn jupyter



Installation of S-timator
~~~~~~~~~~~~~~~~~~~~~~~~~

*S-timator* is on the Python Package Index (pypi), so , after installing the
required libraries, (``Python``, ``numpy``, ``scipy``,
``matplotlib``, ``pandas``, ``seaborn`` and ``pip``) the easiest way to install *S-timator* is
with ``pip``::

    $ pip install stimator

Or, in case `Anaconda/Miniconda` was installed::

    $ conda install -c aeferreira stimator



