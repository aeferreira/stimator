# S-timator

**_ODE dynamical systems modelling in Python_**

## Overview

`S-timator` is a Python module to analyse **ODE-based models**

These models are often found in many scientific fields, particularly in Physics, Chemistry, Biology and
Engineering. They are useful representations of the *dynamics* of systems.

Some of the features that `S-timator` offers are

- [A mini language used to describe models](models.md): models can be input as plain text 
  following a very simple and human-readable language.
- [numerical solution](solving.md) of ODE's, [parameter scanning](solving.md#parameter-scanning).
- [Parameter estimation](par_estimation.md) and *model selection*: given experimental data in
  the form of time series and constrains on model operating ranges,
  built-in numerical optimizers can find parameter values and assist you in the
  experimental design for model selection.

For a brief introduction check out the [basic features](basic_features.md) mini-tutorial.

## Installation

S-timator requires Python >= 3.8.

It is recomended that S-timator is installed in a virtual environment, for example,
using [venv environments](https://docs.python.org/3/library/venv.html) or [conda environments](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html).

The latest stable version of *S-timator* can be installed from the [Python Package Index](https://pypi.python.org/pypi), with `pip`:

``` bash
python -m pip install stimator
```

All the necessary requirements will be installed.

However, `S-timator` takes advantage of the "Scientific Python ecossystem", a set of
Python libraries that brings high-performance scientific computing to the Python
programming language.

The "scientific Python" requirements are:

- Python, version 3.8 or above
- numpy
- scipy
- matplotlib
- sympy
- Jupyter (optional but recommended to run the notebook examples)

Known "scientific distributions" of the Python language are:

- [Anaconda distribution](https://www.anaconda.com/products/distribution)
(or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) followed by the necessary `conda install` commands)
- [Python (x,y)](http://python-xy.github.io/)`
- Enthought distribution through the [Enthought Deployment Manager](https://assets.enthought.com/downloads/)

If using miniconda, the requirements can be installed by running
```
conda install numpy scipy matplotlib sympy jupyter
python -m pip install stimator
```

