# S-timator

**_ODE dynamical systems modelling in Python_**

## Overview

`S-timator` is a Python module to analyse **ODE-based models**

These models are often found in many scientific fields, particularly in Physics, Chemistry, Biology and
Engineering. They are useful representations of the *dynamics* of systems.

Some of the features offered by `S-timator` are

- [A mini language used to describe models](models.md): models can be input as plain text 
  following a very simple and human-readable language.
- [numerical solution](solving.md) of ODE's, [parameter scanning](solving.md#parameter-scanning).
- [Parameter estimation](par_estimation.md) and *model selection*: given experimental data in
  the form of time series and constrains on model operating ranges,
  built-in numerical optimizers can optimize parameter values.

For a brief introduction check out the [basic features](basic_features.md) mini-tutorial.

## Installation

S-timator requires Python >= 3.8.

It is recomended that S-timator is installed in a virtual environment, for example,
using [venv environments](https://docs.python.org/3/library/venv.html) or [conda environments](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html).

The latest stable version of *S-timator* can be installed from the [Python Package Index](https://pypi.python.org/pypi), with `pip`:

``` bash
python -m pip install stimator
```

All the necessary requirements will be installed. These include popular Python modules for scientific computing
such as `numpy`, `scipy`, `matplotlib`, and `sympy`

However, `S-timator` may take advantage of having other packages from the "Scientific Python ecossystem", installed.

For example, `S-timator` may be used form Jupyter notebooks, requiring the installation of a platform to render and run such documents.

As another example, the solutions of ODE systems may be converted to popular _DataFrame_ types used in used in Python, such as Pandas or Polars dataframes, in which case these packages nedd to be installed.

Therefore, it is convenient to install `S-timator` over  "scientific distributions" of the Python language.

Examples are:

- [Anaconda distribution](https://www.anaconda.com/products/distribution)
(or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) followed by the necessary `conda install` commands)
- [Python (x,y)](http://python-xy.github.io/)`
- Enthought distribution through the [Enthought Deployment Manager](https://assets.enthought.com/downloads/)

If using miniconda, `S-timator` can be installed by running
```
conda install numpy scipy matplotlib sympy jupyter
python -m pip install stimator
```

## A basic example

The following example showcases the how a solution of an ODE system can be obtained and plotted:

```py
import stimator as st
st.style.use(['st-seaborn-whitegrid', 'seaborn-talk'])

model="""

"""
```
