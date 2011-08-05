from stimator import *
from time import time, sleep
from numpy import append, linspace

#---------------- EXAMPLE 1
m1 = ui.load_model("glxs_simple.mdl")
print m1
solution1 = solve(m1, tf = 4030.0)

print '--- Last time point ----'
print 'At t =', solution1.t[-1]
for x,value in solution1.last:
    print "%-8s= %f" % (x, value)
#---------------- EXAMPLE 2 
m2 = ui.load_model("branched.mdl")
times = append(linspace(0.0,5.0,500),linspace(5.0,10.0, 500))
solution2 = solve(m2, tf = 10.0, times=times)
#---------------- EXAMPLE 3 
m3 = ui.load_model("calcium.mdl")
solution3 = solve(m3, tf = 8.0, npoints = 2000)
#---------------- EXAMPLE 4
m4 = ui.load_model("rossler.mdl")
solution4 = solve(m4, tf = 100.0, npoints = 2000, outputs="x1 x2 x3")

def transformation(vars,t):
    if t > 40.0:
        return (vars[0]-5.0, vars[1], vars[2])
    else:
        return (-5.0, vars[1], vars[2])

solution4.apply_transf(transformation)

plot ([solution1, solution2, solution3, solution4], figure=ui.new_figure())

#---------------- PARAMETER SCANNING
m = ui.load_model("calcium.mdl")

s = Solutions("CICR model: Effect of stimulus on citosolic calcium")

time0 = time()
print '\n--- Starting parameter scanning...'
for stimulus in 0.0, 0.2, 0.4, 0.78:
    m.B = stimulus
    s += solve(m, tf = 6.0, npoints = 10000, title = 'stimulus = %g'% (m.B), outputs="Ca CaComp")

print 'done in', time()-time0, 's'

plot(s,figure=ui.new_figure(), ynormalize = True)
