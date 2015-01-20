.. _installing:

Installing and getting started
------------------------------


Requirements
~~~~~~~~~~~~

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
~~~~~~~~~~~~

*S-timator* is on the Python Package Index (pypi), so , after installing the
required libraries, (``Python``, ``numpy``, ``scipy``,
``matplotlib`` and ``pip``) the easiest way to install *S-timator* is
with ``pip``::

    $ pip install stimator



