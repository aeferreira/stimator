from stimator import read_model, readTCs, solve
from stimator.deode import DeODESolver

mdl = """# Example file for S-timator
title Example 2

vin  : -> x1     , rate = k1
v2   : x1 ->  x2 , rate = k2 * x1
vout : x2 ->     , rate = k3 * x2
k1 = 1
k2 = 2
k3 = 1
init = state(x1=0, x2=0)
!! x2
find k1  in [0, 2]
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

solver = DeODESolver(m1,optSettings, timecourses)
solver.Solve()
print solver.reportResults()
solver.draw()

m2 = m1.copy()
best = solver.optimum.parameters
best = [(n,v) for n,v,e in best]
m2.update(best)
solve(m2, tf=20.0).plot(show=True)
