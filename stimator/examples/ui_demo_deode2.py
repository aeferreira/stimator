"""S-timator : DEMO of deode module."""

from stimator import *
from stimator.deode import DeODESolver

print __doc__
print
print """The deode module combines ODE solving with DE (differential evolution)
"""

m1 = ui.load_model("glxs_hta.mdl")
optSettings = m1['optSettings']
tc = readTCs(m1, ".", verbose=True)

solver = DeODESolver(m1,optSettings, tc)
#solver.Solve()
while solver.exitCode ==0:
    ui.checkpoint()
    solver.computeGeneration()
solver.finalize()
print solver.reportResults()

fig1 = ui.new_figure()
solver.draw(fig1)

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

#solver.Solve()
while solver.exitCode ==0:
    ui.checkpoint()
    solver.computeGeneration()
solver.finalize()
print solver.reportResults()

fig2 = ui.new_figure()
solver.draw(fig2)
