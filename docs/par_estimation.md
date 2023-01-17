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

# Parameter estimation.

The **estimation.py** module combines ODE solving with the DE (differential evolution) genetic optimizer.

```{code-cell} ipython3
import stimator as st
```

##  Linear pathway with three reactions

```{code-cell} ipython3

# ----------- Model ------------------------

mdl = """# Example file for S-timator
title Example 1

vin  : -> x1     , rate = k1
v2   : x1 ->  x2 , rate = k2 * x1
vout : x2 ->     , rate = k3 * x2

init : x1=0, x2=0

find k1 in [0, 2]
find k2 in [0, 2]
find k3 in [0, 2]

!! x1 x2

popsize = 60     # population size in GA
"""

m1 = st.read_model(mdl)

# ----------- Time course -------------------

example_data = """
t   x1   x2
0   0   0
2   1.403812093   0.48351624
4   1.528870297   1.483289613
6   1.917963699   2.039584833
8   2.028998372   2.826410056
10   1.978326655   3.106415222
12   2.143692636   3.060669986
14   2.289572191   3.231815374
16   2.019850835   3.310127564
18   1.977904321   3.098886165
20   2.126776717   3.463202683
"""
```

```{code-cell} ipython3
best = m1.estimate(timecourses=example_data)

print(best)
```

One can update the model parameters to the best fit values and obtain the same timecourse

```{code-cell} ipython3
m2 = m1.copy()
bestpars = [(n,v) for n,v,e in best.parameters]
m2.setp(bestpars)
dict(bestpars)
```

A bit of styling of the plots:

```{code-cell} ipython3
%matplotlib inline
from matplotlib import pyplot as plt
st.style.use(['st-seaborn-whitegrid', 'seaborn-talk'])
```

```{code-cell} ipython3
# plotting side by side
_, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

best.plot(ax=ax1, palette='Dark2', xlabel='time')
m2.solve(tf=20.0).plot(ax=ax2, palette='Dark2', xlabel='time')

plt.show()
```

##  Glyoxalase system

+++

### An example with **two time courses**

```{code-cell} ipython3
mdl = """
title example 2: Glyoxalase system in L. Infantum

glx1 : HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
#glx1 : HTA -> SDLTSH, V*HTA/(Km1 + HTA), V=2.57594e-05
glx2 : SDLTSH ->,     V2*SDLTSH/(Km2 + SDLTSH)

#find glx1.V  in [0.00001, 0.0001]
find V1  in [0.00001, 0.0001]

Km1 = 0.252531
find Km1 in [0.01, 1]

V2  = 2.23416e-05
find V2 in [0.00001, 0.0001]

Km2 = 0.0980973
find Km2 in (0.01, 1)

init : (SDLTSH = 7.69231E-05, HTA = 0.1357)

timecourse TSH2a.txt
timecourse TSH2b.txt
"""
m1 = st.read_model(mdl)
print(mdl)
```

```{code-cell} ipython3
tcdir = st.get_examples_path()

optimum = m1.estimate(tc_dir=tcdir, names=['SDLTSH', 'HTA'])
print(optimum)
```

```{code-cell} ipython3
#plt.style.use('bmh')

_, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey='row')
optimum.plot(0, ax=ax1, ylabel='conc (microM)')
optimum.plot(1, ax=ax2, xlabel='t (s)');
```

-----------

### An example with an *unknown initial value*

```{code-cell} ipython3
m2 = m1.copy()

# Assume init.HTA is uncertain
m2.init.HTA.set_bounds((0.05,0.25))

# do not estimate Km1 and Km2, just to help the analysis
m2.parameters.Km1.set_bounds(None)
m2.parameters.Km2.set_bounds(None)
m2.parameters.Km1 = 0.252531
m2.parameters.Km2 = 0.0980973


# VERY IMPORTANT:
# only one time course can be used: 
# cannot fit one initial value using several timecourses!

best = m2.estimate('TSH2a.txt',
                  names=['SDLTSH', 'HTA'], tc_dir=tcdir,
                  opt_settings=dict(pop_size=60))

print(best)
```

```{code-cell} ipython3
best.plot();
```

### An example with a transformation

```{code-cell} ipython3
mtransf = st.read_model("""
title example 2, fitting a transformation

glx1 : HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
#glx1 : HTA -> SDLTSH, V*HTA/(Km1 + HTA), V=2.57594e-05
glx2 : SDLTSH ->,     V2*SDLTSH/(Km2 + SDLTSH)

#find glx1.V  in [0.00001, 0.0001]
find V1  in [0.00001, 0.0001]

Km1 = 0.252531
find Km1 in [0.01, 1]

V2  = 2.23416e-05
find V2 in [0.00001, 0.0001]

Km2 = 0.0980973
find Km2 in (0.01, 1)

~sdlx2 = 2 * SDLTSH # the transformation to fit

!! sdlx2

init : (SDLTSH = 7.69231E-05, HTA = 0.1357)

""")

optimum = mtransf.estimate(tc_dir=tcdir,
                           timecourses='tc_double.txt',
                           names=['sdlx2', 'SDLTSH', 'HTA'])

print(optimum)
optimum.plot()
plt.show()
```
