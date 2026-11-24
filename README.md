# S-timator

Analysis of ODE models with focus on model selection and parameter estimation.

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

## Requirements

*S-timator* supports Python versions 3.12+.

*S-timator* depends on the "scientific python stack". The **mandatory**
requirements for *S-timator* are the following libraries:

- `Python` 3.12 or higher
- `matplotlib` 2.0 or higher
- `numpy`
- `scipy`
- `sympy`
- `packaging`
- `dotmap`


## Installation

with `pip`

```
python -m pip install stimator
```

The installation of these Python libraries is optional, but strongly recommended:

- `pandas` or `polars`, to treat time courses as DataFrames
- `Jupyter` or `marimo` notebook platforms

## Short example

```python
print("placeholder")
```
