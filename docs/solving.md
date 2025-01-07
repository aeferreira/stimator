# Solving ODE models.

This notebook shows how to use 4 of the most common **S-timator** functions:

- `read_model()`, reads a _string_ that conforms to a model description language, returning a `Model` object
- `solve()`, computes a solution of the ODE system associated with a model.
- `scan()`, calls `Model.solve()` several times, scanning a model parameter in a range of values.
- `plot()`, draws a graph of the results returned from `solve()` or `scan()`.

```py exec="true" session="solving" source="above"
import stimator as st
from stimator import Solutions
from stimator.plots import prepare_grid
import stimator.examples.models as models
```

Before we begin, a bit of styling of the plots:

```py exec="true" session="solving" source="above"
from io import StringIO # markdown-exec: hide
from matplotlib import pyplot as plt
st.style.use(['st-seaborn-whitegrid', 'seaborn-talk'])
```

## Example 1: Glyoxalase system

```py exec="true" session="solving" source="above" result="txt"
mdl = models.glyoxalases.text
print(mdl)
m1 = st.read_model(mdl)
```

```py exec="true" session="solving" source="above" html="1"
f, ax = plt.subplots() # markdown-exec: hide
s = m1.solve(tf=4030.0)

s.plot(xlabel='$t$ (min)');

buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

```py exec="true" session="solving" source="above" result="txt"
print(f'Last time point at t = {s.t[-1]}')
for name, value in s.last.items():
    print(f"{name:8s} = {value:.3f}")
```

## Example 2: Branched pathway

```py exec="true" session="solving" source="above" result="txt"
mdl = models.branched.text

print(mdl)
```

```py exec="true" session="solving" source="above" html="1"
plt.close() # markdown-exec: hide
f, ax = plt.subplots() # markdown-exec: hide
m2 = st.read_model(mdl)

from numpy import append, linspace
times = append(linspace(0.0, 5.0, 500), linspace(5.0, 10.0, 500))

m2.solve(tf=10.0, times=times).plot();
buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

## Example 3: Calcium spikes: CICR model

```py exec="true" session="solving" source="above" result="txt"
mdl = models.ca.text

print(mdl)
```

```py exec="true" session="solving" source="above" html="1"
plt.close() # markdown-exec: hide
f, ax = plt.subplots() # markdown-exec: hide
#chaining functions...
st.read_model(mdl).solve(tf=8.0, npoints=2000).plot();

buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

## Example 4: Rossler chaotic system

```py exec="true" session="solving" source="above" result="txt"
mdl = models.rossler.text
print (mdl)
```

``` py exec="true" session="solving" source="above" html="1"
plt.close() # markdown-exec: hide
f, ax = plt.subplots() # markdown-exec: hide
m4 = st.read_model(mdl)

s = m4.solve(tf=100.0, npoints=2000, outputs="x1 x2 x3".split())

def transformation(vars, t):
    if t > 40.0:
        return (vars[0] - 5.0, vars[1], vars[2])
    else:
        return (-5.0, vars[1], vars[2])

s.apply_transf(transformation)

s.plot();

buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

## Example 5: Lorentz system (sensitivity to initial conditions)

```py exec="true" session="solving" source="above" result="txt"
mdl = models.lorentz.text
print (mdl)
```

``` py exec="true" session="solving" source="above" html="1"
plt.close() # markdown-exec: hide
m5 = st.read_model(mdl)

ivs = {'init.x':(1.0, 1.01, 1.02)}
titles = [f'$x(0)$ = {iv}' for iv in ivs['init.x']]

s = m5.scan(ivs, tf=25.0, npoints=20000, outputs=['x'], titles=titles)

f, ax = plt.subplots() # markdown-exec: hide
s.one_plot(what='x', label_fmt='$title', title=m5.metadata['title']);

buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

## Parameter scanning

### Example 6: parameter scanning in the CICR model

```py exec="true" session="solving" source="above" html="1"
plt.close() # markdown-exec: hide
m = st.read_model("""
title Calcium Spikes
v0         = -> Ca, 1
v1         = -> Ca, k1*B*step(t, 1.0), k1 = 7.3

B          = 0.4
export     = Ca ->  , 10 ..
leak       = CaComp -> Ca, 1 ..

!! Ca

v2         = Ca -> CaComp, 65 * Ca**2 / (1+Ca**2)
v3         = CaComp -> Ca, 500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)

init       : (Ca = 0.1, CaComp = 0.63655)
""")

bvalues = (0.0, 0.1, 0.2, 0.25, 0.28, 0.29, 0.3,
           0.35, 0.4, 0.45, 0.5, 0.6, 0.75, 0.8, 0.9, 1.0)
titles = [f'$\\beta$ = {b:g}' for b in bvalues]

s = m.scan({'B': bvalues}, tf=8.0, npoints=1000, titles=titles)
suptitlegend="Dynamics of cytosolic $Ca^{2+}$ as a function of stimulus"

f, axs = prepare_grid(s, figsize=(16,16), constrained_layout=True)

s.plot(ylim=(0,1.5), axs=axs, legend=False, xlabel='$t$ (min)')

f.suptitle(suptitlegend, fontsize=20);

buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```

Several time courses in the same plot

```py exec="true" session="solving" source="above" html="1"
plt.close() # markdown-exec: hide
f, ax = plt.subplots() # markdown-exec: hide
sols = Solutions([s[i] for i in range(0, len(s), 3)])
sols.one_plot(ylim=(0,1.5),
              legend='out',
              title='CICR model of Calcium Spikes',
              xlabel='$t$ (min)',
              label_fmt='$title',
              palette='tab20');

buffer = StringIO() # markdown-exec: hide
f.savefig(buffer, format="svg") # markdown-exec: hide
print(buffer.getvalue()) # markdown-exec: hide
```
