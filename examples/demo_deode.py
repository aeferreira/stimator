import sys
import os.path

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir))

from time import time
from stimatordev import *
from stimatordev.deode import DeODESolver

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

m1 = Model("Glyoxalase system in L.infantum")
m1.glo1 = react("HTA -> SDLTSH", rate = "V1*HTA/(Km1 + HTA)")
m1.glo2 = react("SDLTSH -> "   , rate = "V2*SDLTSH/(Km2 + SDLTSH)")
m1.V1  = 2.57594e-05
m1.Km1 = 0.252531
m1.V2  = 2.23416e-05
m1.Km2 = 0.0980973
m1.V1.uncertainty(0.00001, 0.0001)
m1.Km1.uncertainty(0.01, 1)
m1.V2.uncertainty(0.00001, 0.0001)
m1.Km2.uncertainty(0.01, 1)
m1.init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)
#print m1

optSettings={'genomesize':80, 'generations':200}
timecourses = readTimeCourses(['TSH2a.txt', 'TSH2b.txt'], '.', (0,2,1))

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
timecourses = readTimeCourses(['TSH2a.txt'], '.', (0,2,1))

solver = DeODESolver(m2,optSettings, timecourses)

time0 = time()

solver.Solve()

print "Optimization took %f s"% (time()-time0)

print
print '---------------------------------------------------------'
print "Results for %s\n" % m2.title
print reportResults(solver)
