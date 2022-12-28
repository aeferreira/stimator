# S-timator: ODE dynamical systems modelling in python

`S-timator` is a Python module to analyse **ODE-based models**
(also known as *dynamic* or *kinetic* models).

These models are often found in many scientific fields, particularly in Physics, Chemistry, Biology and
Engineering.

Some of the features that `S-timator` offers are

- **A mini language used to describe models**: models can be input as plain text 
  following a very simple and human-readable language.
- **Basic analysis**: numerical solution of ODE's, parameter scanning.
- **Parameter estimation** and **model selection**: given experimental data in
  the form of time series and constrains on model operating ranges,
  built-in numerical optimizers can find parameter values and assist you in the
  experimental design for model selection.

`S-timator` is in an **alpha stage**: many new features will be available soon.

For a brief introduction to the package, you can read the :ref:`Basic use <introduction>`.

Much more detail can be found in the :ref:`tutorial <tutorial>`.

You can also browse the :ref:`API reference <api_ref>` to see the kind of tools that are available.

## Installation: Python and "Scientific" Python

`S-timator` is a (pure) Python package that is installable
from the [Python Package Index](https://pypi.python.org/pypi).

However, `S-timator` relies heavily on the "Scientific Python ecossystem", a set of
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


### Requirements

`S-timator` supports Python versions 3.6 and above.

*S-timator* depends on a "scientific python stack". The **mandatory**
requirements for *S-timator* are the following libraries:

- `numpy`
- `scipy`
- `matplotlib`
- `sympy`


Other Python libraries which are optional, but strongly recommended:

- The *Jupyter* package: some *S-timator* examples are provided
  as Jupyter notebooks.

The `Anaconda Python Distribution, from "Continuum Analytics" <https://store.continuum.io/cshop/anaconda/>`_
is a convenient Python distribution that includes all of the above requirements.

If the full `Anaconda Distribution` is too heavy on disk space, from the same company, the `Miniconda <http://conda.pydata.org/miniconda.html>`_ "slim" distribution
is also an alternative. In this case, run the necessary ``conda install``'s of the requirements, after installing Miniconda:

```
conda install numpy scipy matplotlib sympy jupyter
```

### Installation of S-timator

*S-timator* is on the Python Package Index (pypi), so , after installing the
required libraries, (``Python``, ``numpy``, ``scipy``,``matplotlib`` and ``pip``)
the easiest way to install *S-timator* is with `pip`:`

```
pip install stimator
```










