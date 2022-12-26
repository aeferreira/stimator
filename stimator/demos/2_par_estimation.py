from stimator import read_model, __version__, Solution
from matplotlib import pyplot as plt

print('S-timator version', __version__)

# --- example 1 --------------------

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
tc = Solution.read_str(example_data)

mdl = """# Example model
title Model for parameter estimation

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

model = read_model(mdl)

best = model.estimate(timecourses=tc)

print(best)

print('--- Modifying model ---')
m2 = model.copy()
bestpars = [(n, v) for n, v, _ in best.parameters]
print(dict(bestpars))

m2.setp(bestpars)

# plot side by side
_, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

best.plot(ax=ax1, palette='Dark2', xlabel='t')
m2.solve(tf=20.0).plot(ax=ax2, palette='Dark2', xlabel='t')

plt.show()
