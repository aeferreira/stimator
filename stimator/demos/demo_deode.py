"""S-timator : DEMO of deode module."""

from time import time
from stimator import *
from stimator.deode import DeODESolver

print __doc__
print
print """The deode module combines ODE solving with DE (differential evolution)
"""

def reportResults(solver):
    reportText = ""
    sections = [solver.optimum[s] for s in ['parameters', 'optimization', 'timecourses']]
    for section in sections:
        reportText += "--- %-20s -----------------------------\n" % section['name'].upper()
        if section['header']:
            reportText += '\t\t'.join(section['header'])+'\n'
        reportText += "\n".join([section['format'] % i for i in section['data']])
        reportText += '\n\n'
    return reportText

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
timecourses = readTimeCourses(['TSH2a.txt', 'TSH2b.txt'], '../../models', (0,2,1))

solver = DeODESolver(m1,optSettings, timecourses)

time0 = time()

solver.Solve()

print "Optimization took %f s"% (time()-time0)

print
print '---------------------------------------------------------'
print "Results for %s\n" % m1.title
print reportResults(solver)

#--- an example with unknown initial values --------------------

m2 = m1.clone()

# Now, assume init.HTA is uncertain
m2.init.HTA.uncertainty(0.05,0.25)
# do not estimate Km1 and Km2 to help the analysis
m2.Km1.uncertainty(None)
m2.Km2.uncertainty(None)

optSettings={'genomesize':60, 'generations':200}

## VERY IMPORTANT:
## only one time course can be used: 
## cannot fit one uncertain initial value to several timecourses!!!
timecourses = readTimeCourses(['TSH2a.txt'], '../../models', (0,2,1))

solver = DeODESolver(m2,optSettings, timecourses)

time0 = time()

solver.Solve()

print "Optimization took %f s"% (time()-time0)

print
print '---------------------------------------------------------'
print "Results for %s\n" % m2.title
print reportResults(solver)
