---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Solving ODE models (`dynamics.py`).

```{code-cell} ipython3
%matplotlib inline
from matplotlib import pyplot as plt
plt.style.use(['seaborn-whitegrid',
              {'xaxis.labellocation': 'right',
               'legend.frameon': True,
               'legend.facecolor': 'white'}])
```

+++

This notebook shows how to use 4 of the most common **S-timator** functions:

- `read_model()`, reads a _string_ that conforms to a model description language, returning a `Model` object
- `solve()`, computes a solution of the ODE system associated with a model.
- `scan()`, calls Model.solve() several times, scanning a model parameter in a range of values.
- `plot()`, draws a graph of the results returned from `solve()` or `scan()`.

```{code-cell} ipython3
from stimator import read_model, examples, Solutions
from stimator.plots import prepare_grid
```

## Example 1: Glyoxalase system

```{code-cell} ipython3
mdl = examples.models.glyoxalases.text
print(mdl)
m1 = read_model(mdl)

s = m1.solve(tf=4030.0)
s.plot()

print('==== Last time point ====')
print('At t = %g'% s.t[-1])
for x in s.last:
    print("%-8s= %f" % (x, s.last[x]))
```

## Example 2: Branched pathway

```{code-cell} ipython3
from numpy import append, linspace
mdl = examples.models.branched.text

print(mdl)

m2 = read_model(mdl)

times = append(linspace(0.0, 5.0, 500), linspace(5.0, 10.0, 500))

m2.solve(tf=10.0, times=times).plot()
```

## Example 3: Calcium spikes: CICR model

```{code-cell} ipython3
mdl = examples.models.ca.text

print(mdl)

#chaining functions...
read_model(mdl).solve(tf=8.0, npoints=2000).plot()
```

## Example 4: Rossler chaotic system

```{code-cell} ipython3
mdl = examples.models.rossler.text; print (mdl)
m4 = read_model(mdl)

s = m4.solve(tf=100.0, npoints=2000, outputs="x1 x2 x3".split())

def transformation(vars, t):
    if t > 40.0:
        return (vars[0] - 5.0, vars[1], vars[2])
    else:
        return (-5.0, vars[1], vars[2])

s.apply_transf(transformation)

s.plot();
```

## Example 5: Lorentz system (sensitivity to initial conditions)

```{code-cell} ipython3
mdl = examples.models.lorentz.text
print (mdl)
m5 = read_model(mdl)

ivs = {'init.x':(1.0, 1.01, 1.02)}
titles = ['$x(0)$ = %g' % v for v in ivs['init.x']]

s = m5.scan(ivs, tf=25.0, npoints=20000, outputs=['x'], titles=titles)

f, ax = plt.subplots(figsize=(12, 8))
s.one_plot(what='x', ax=ax, label_fmt='$title');
f.suptitlegend=m5.metadata['title'],
```

## Example 6: parameter scanning in the CICR model

```{code-cell} ipython3
m = read_model("""
title Calcium Spikes
v0         = -> Ca, 1
v1         = -> Ca, k1*B*step(t, 1.0)
k1         = 7.3
B          = 0.4
export     = Ca ->  , 10 ..
leak       = CaComp -> Ca, 1 ..
!! Ca
v2         = Ca -> CaComp, 65 * Ca**2 / (1+Ca**2)
v3         = CaComp -> Ca, 500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)
init       : (Ca = 0.1, CaComp = 0.63655)""")

bvalues = (0.0, 0.1, 0.2, 0.25, 0.28, 0.29, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.75, 0.8, 0.9, 1.0)
titles = [f'$\\beta$ = {b:g}' for b in bvalues]

s = m.scan({'B': bvalues}, tf=8.0, npoints=1000, titles=titles)
suptitlegend="Dynamics of cytosolic $Ca^{2+}$ as a function of stimulus"

f, axs = prepare_grid(s, figsize=(16,16))

s.plot(ylim=(0,1.5), axs=axs, legend=False, xlabel='$t$ (min)')
f.suptitlegend = suptitlegend
```

Several time courses in the same plot

```{code-cell} ipython3
sols = Solutions([s[i] for i in range(0, len(s), 3)])
sols.one_plot(ylim=(0,1.5),
              legend='out',
              xlabel='$t$ (min)',
              label_fmt='$title',
              palette='tab20');
```
