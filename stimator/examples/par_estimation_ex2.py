from stimator import read_model, readTCs
from stimator.deode import DeODEOptimizer

mdl = """# Example file for S-timator
title Example 2

vin  : -> x1     , rate = k1
v2   : x1 ->  x2 , rate = k2 * x1
vout : x2 ->     , rate = k3 * x2

init : x1=0, x2=0

find k1 in [0, 2]
find k2 in [0, 2]
find k3 in [0, 2]

timecourse ex2data.txt
generations = 200   # maximum generations for GA
genomesize = 60     # population size in GA
"""
m1 = read_model(mdl)
print mdl

optSettings={'genomesize':60, 'generations':200}
timecourses = readTCs(['ex2data.txt'], verbose=True)

optimizer = DeODEOptimizer(m1,optSettings, timecourses)
optimizer.run()
best = optimizer.optimum
print best.info()
best.plot()

m2 = m1.copy()
bestpars = [(n,v) for n,v,e in best.parameters]
m2.update(bestpars)
#solve(m2, tf=20.0).plot(show=True)
m2.solve(tf=20.0, outputs='x2').plot(show=True)
