"""S-timator : DEMO of deode module."""

from stimator import *
from stimator.deode import DeODESolver
import pylab as pl

print __doc__
print
print """The deode module combines ODE solving with DE (differential evolution)
"""

m1 = read_model("""
title Glyoxalase system in L. Infantum

#reactions (with stoichiometry and rate)
glx1 : HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
glx2 : SDLTSH ->,     V2*SDLTSH/(Km2 + SDLTSH)

find V1  in [0.00001, 0.0001]
find Km1 in [0.01, 1]
find V2  in [0.00001, 0.0001]
find Km2 in [0.01, 1]
init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)
""")

#~ print m1

optSettings={'genomesize':80, 'generations':200}
timecourses = readTCs(['TSH2a.txt', 'TSH2b.txt'], names = ['SDLTSH', 'HTA'], verbose=True)

solver = DeODESolver(m1,optSettings, timecourses)
solver.Solve()

print solver.reportResults()

fig1 = pl.figure()
solver.draw(fig1)

redsols = solver.optimum.optimum_tcs
redsols.saveTimeCoursesTo(['TSH2a_pred.txt', 'TSH2b_pred.txt'], verbose=True)
## for tcs in solver.fitted_tcs:
##     plot(tcs, superimpose=True)


#--- an example with unknown initial values --------------------

m2 = m1.clone()

# Now, assume init.HTA is uncertain
m2.init.HTA.uncertainty(0.05,0.25)
# do not estimate Km1 and Km2 to help the analysis
m2.Km1.uncertainty(None)
m2.Km2.uncertainty(None)
m2.Km1 = 0.252531
m2.Km2 = 0.0980973

optSettings={'genomesize':60, 'generations':200}

## VERY IMPORTANT:
## only one time course can be used: 
## cannot fit one uncertain initial value to several timecourses!!!
timecourses = readTCs(['TSH2a.txt'], '.', names = ['SDLTSH', 'HTA'], verbose=True)

solver = DeODESolver(m2,optSettings, timecourses)
solver.Solve()

print solver.reportResults()

fig2 = pl.figure()
solver.draw(fig2)

pl.show()
