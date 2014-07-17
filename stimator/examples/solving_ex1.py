from stimator import read_model, solve
import pylab as pl

mdl = """# Example file for S-timator
title Example 1

#reactions (with stoichiometry and rate)
vin  : -> x1     , rate = k1
v2   : x1 ->  x2 , rate = k2 * x1
vout : x2 ->     , rate = k3 * x2

#parameters and initial state
k1 = 1
k2 = 2
k3 = 1
init = state(x1=0, x2=0)

#filter what you want to plot
!! x1 x2"""

m = read_model(mdl)

print '========= model ========================================'
print mdl
print '--------------------------------------------------------'

s1 = solve(m, tf=5.0)
s1.plot()

pl.show()