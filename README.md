# S-timator

**_ODE dynamical systems modelling in Python_**

## Overview

`S-timator` is a Python module to analyse **SODE-based models**

These models are often found in many scientific fields, particularly in Physics, Chemistry, Biology and Engineering.  
They are useful representations of the *dynamics* of systems.

Some of the features offered by `S-timator` are

- A **mini language** used to describe models: models can be input as plain text 
  following a very simple and human-readable language.
- **Numerical solution** of SODE and **parameter scanning**.
- **Parameter estimation** and *model selection*: given experimental data in  
the form of time series and constrains on model operating ranges,
built-in numerical optimizers can optimize parameter values.

For a brief introduction check out the mini-tutorial in the documentation.

## Installation

S-timator requires Python >= 3.12.

S-timator (module `stimator`) is a [package on Python Package Index](https://pypi.org/project/stimator/). As such, it can be installed in a working environment by using one of the following package managers:

- [pip](https://pip.pypa.io/en/stable/)  
`python -m pip install stimator`, preferably in a `venv` or `conda` environment.
- [uv](https://docs.astral.sh/uv/)  
`uv add stimator` to add it to a uv managed project or `uv pip install stimator` to add it to any virtual environment
- [pixi](https://pixi.prefix.dev/latest/)  
`pixi add --pypi stimator` to add it to a pixi managed project or workspace

## Requirements

All the necessary requirements will be installed. These are

- `numpy`
- `scipy`
- `matplotlib`
- `sympy`
- `packaging`
- `dotmap`

It is recommended to pre-install these requirements as conda packages, by using [pixi](https://pixi.prefix.dev/latest/) or the [Anaconda/Miniconda distribution](https://www.anaconda.com/download/success?reg=skipped-miniconda&path=miniconda)

`S-timator` may take advantage of having other tools from the "Scientific Python ecossystem". For example, `S-timator` may be used form [Jupyter](https://jupyter.org/) or [Marimo](https://marimo.io/) notebooks.

As another example, the solutions of ODE systems may be converted to popular _DataFrame_ types, such as Pandas or Polars dataframes.

## A warm-up example

The following example showcases the how the solution of an ODE system can be obtained:

```python
import pandas as pd
import stimator as st

model = st.read_model("""
title A simple open chemical system

vin: -> A, rate = 0.5
r1: A -> B, rate = k1 * A
r2: B -> C, rate = k2 * B - k3 * C
vout: C ->, rate = 0.8 * C

k1 = 0.1
k2 = 2
k3 = 1

init: A = 0, B = 0, C = 0
""")

# compute the solution of the SODE up to time 80
solution = model.solve(tf=80)

# show the solution as a table
print(pd.DataFrame(solution.as_dict()).set_index("t"))

# Report the last point of the solution
print(f"""
Steady state:
{solution.last}
""")
```

```
                  A         B         C
t                                      
0.000000   0.000000  0.000000  0.000000
0.160321   0.079521  0.000578  0.000059
0.320641   0.157778  0.002106  0.000409
0.480962   0.234789  0.004358  0.001207
0.641283   0.310576  0.007186  0.002519
...             ...       ...       ...
79.358717  4.998212  0.562253  0.624709
79.519038  4.998240  0.562257  0.624714
79.679359  4.998268  0.562261  0.624718
79.839679  4.998296  0.562264  0.624723
80.000000  4.998323  0.562268  0.624727

[500 rows x 3 columns]

Steady state:
{'A': 4.998322658118175, 'B': 0.5622681719167945, 'C': 0.6247272619351981}
```

The solution could easily be plotted with

```python
solution.plot(xlabel="time")
```