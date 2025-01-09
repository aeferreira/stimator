# Parameter estimation.

Given experimental observations, *parameter estimation* is the tuning of the parameters of an ODE model so that its solutions fit the observations as best as possible, usually using a least-squares measure for fitting.

The **estimation.py** module combines ODE solving with the DE (differential evolution) genetic optimizer.

As indicated in the [basic features](basic_features.md) page, we start by importing `stimator`, and tweeking a bit the matplotlib style for plotting:

```py exec="true" source="above" session="parst"
import stimator as st
from io import StringIO # markdown-exec: hide
from matplotlib import pyplot as plt
st.style.use(['st-seaborn-whitegrid', 'seaborn-talk'])
```

##  A simple example

```py exec="true" source="above" session="parst"

mdl = """title Example 1

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

```py exec="true" source="above" session="parst" result="txt"
best = m1.estimate(timecourses=example_data)

print(best.progress_report) # markdown-exec: hide
print(best)
```

One can update the model parameters to the best fit values and obtain the same timecourse

```py exec="true" source="above" session="parst" result="txt"
m2 = m1.copy()
bestpars = [(n,v) for n,v,e in best.parameters]
m2.setp(bestpars)
for (pname, pvalue) in bestpars:
    print(f"{pname} = {pvalue}")
```

```py exec="true" source="above" session="parst" html="1"
# plotting side by side
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

best.plot(ax=ax1, palette='Dark2', xlabel='time')
m2.solve(tf=20.0).plot(ax=ax2, palette='Dark2', xlabel='time')
buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

## An example with **two time courses**

### Glyoxalase system

```py exec="true" source="above" session="parst" result="txt"
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

tcdir = st.get_examples_path()

optimum = m1.estimate(tc_dir=tcdir, names=['SDLTSH', 'HTA'])
print(optimum.progress_report) # markdown-exec: hide
print(optimum)
```

```py exec="true" source="above" session="parst" html="1"
plt.close() # markdown-exec: hide

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey='row')
optimum.plot(0, ax=ax1, ylabel='conc (microM)')
optimum.plot(1, ax=ax2, xlabel='t (s)');
buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

### An example with an *unknown initial value*

```py exec="true" source="above" session="parst" result="txt"
m2 = m1.copy()

# Assume init.HTA is uncertain
m2.init.HTA.set_bounds((0.05,0.25))

# do not estimate Km1 and Km2, just to help the analysis
m2.parameters.Km1.set_bounds(None)
m2.parameters.Km2.set_bounds(None)
m2.parameters.Km1 = 0.252531
m2.parameters.Km2 = 0.0980973


# IMPORTANT:
# only one time course can be used: 
# cannot fit one initial value using several timecourses!

best = m2.estimate('TSH2a.txt',
                  names=['SDLTSH', 'HTA'], tc_dir=tcdir,
                  opt_settings=dict(pop_size=60))

print(best.progress_report) # markdown-exec: hide
print(best)
```

```py exec="true" source="above" session="parst" html="1"
plt.close() # markdown-exec: hide
f, ax = plt.subplots() # markdown-exec: hide

best.plot();
buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

### An example with a transformation

```py exec="true" source="above" session="parst" result="txt"
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

print(best.progress_report) # markdown-exec: hide
print(optimum)
```

```py exec="true" source="above" session="parst" html="1"
plt.close() # markdown-exec: hide
f, ax = plt.subplots() # markdown-exec: hide

optimum.plot()
buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```
