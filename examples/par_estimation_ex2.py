from stimator import read_model
from stimator.versions import version_info

print ('Version information')
print (version_info())

mdl = """# Example file for S-timator
title Example 2

vin  : -> x1     , rate = k1
v2   : x1 ->  x2 , rate = k2 * x1
vout : x2 ->     , rate = k3 * x2

init : x1=0, x2=0

find k1 in [0, 2]
find k2 in [0, 2]
find k3 in [0, 2]

!! x2

timecourse ex2data.txt
popsize = 60     # population size in GA
"""
m1 = read_model(mdl)
print mdl

best = m1.estimate()

print best.info()
best.plot()

m2 = m1.copy()
bestpars = [(n,v) for n,v,e in best.parameters]
m2.update(bestpars)

m2.solve(tf=20.0).plot(show=True)
