# S-timator: ODE dynamical systems modelling in python

## Global description

`S-timator` is a Python module to analyse **ODE-based models**
(also known as *dynamic* or *kinetic* models).

These models are often found in many scientific fields, particularly in Physics, Chemistry, Biology and
Engineering.

Some of the features that `S-timator` offers are

- [A mini language used to describe models](models): models can be input as plain text 
  following a very simple and human-readable language.
- **Basic analysis**: [numerical solution](solving_models) of ODE's, [parameter scanning](scanning).
- [Parameter estimation](par_estimation) and **model selection**: given experimental data in
  the form of time series and constrains on model operating ranges,
  built-in numerical optimizers can find parameter values and assist you in the
  experimental design for model selection.

`S-timator` is in an **alpha stage**: many new features will be available soon.

For a brief introduction check out the [basic features](basic-features) mini-tutorial.

## Installation

The latest stable version of *S-timator* can be installed from the [Python Package Index](https://pypi.python.org/pypi), with `pip`:`

```
pip install stimator
```

However, `S-timator` relies heavily on the "Scientific Python ecossystem", a set of
Python libraries that brings high-performance scientific computing to the Python
programming language.

The "scientific Python" requirements are:

- Python, version 3.6 and above
- numpy
- scipy
- matplotlib
- sympy
- Jupyter (optional but recommended to run some of the notebook examples)

In practical terms, this means that the installation of a standard Python
distribution (3.6 or higher) will not be enough to use S-timator. Instead, either the
mandatory dependencies are installed one by one, or, more conveniently, Python
is installed through a "scientific distribution".

One of the following "scientific python" distributions is recommended,
**as they all provide an easy installation of all requirements**:

- [Anaconda distribution](https://www.anaconda.com/products/distribution)
(or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) followed by the necessary `conda install` commands)
- [Python (x,y)](http://python-xy.github.io/)`
- Enthought distribution through the [Enthought Deployment Manager](https://assets.enthought.com/downloads/)

If using miniconda, the requirements can be installed by running
```
conda install numpy scipy matplotlib sympy jupyter
pip install stimator
```

Notice that stimator can be installed in isolated `venv` or `conda` environments.
