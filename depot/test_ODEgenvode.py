import modelparser
#import stimator_timecourse
from numpy import *
import pylab as p
import time

modelText = """
variables: SDLTSH HTA

glx1 : HTA -> SDLTSH, rate = V1*HTA/(Km1 + HTA)
glx2 : SDLTSH -> , rate = V2*SDLTSH/(Km2 + SDLTSH)

find V1 in [0.00001, 0.0001]
find Km1 in [0.01, 1]
find V2 in [0.00001, 0.0001]
find Km2 in [0.01, 1]

timecourse TSH2a.txt
timecourse TSH2b.txt

generations = 200
genomesize = 80

"""


#parameters from optimization using TSH2a.txt and TSH2b.txt
par_names = ('V1', 'Km1', 'V2', 'Km2')
var_names = ('SDLTSH', 'HTA')

V1 = 2.57594e-05
Km1 = 0.252531
V2 = 2.23416e-05
Km2 = 0.0980973
m_Parameters = (V1, Km1, V2, Km2)

header, data = modelparser.readTimeCourseFromFile('../examples/TSH2a.txt')

y0 = copy(data[0, 1:])
t  = data[:, 0]
scale = 10.0*(t[-1]-t[0])
t0 = float(t[0])
times = (t-t0)/scale+t0

parser = modelparser.StimatorParser()
textlines = modelText.split("\n")

parser.parse(textlines)
#modelparser.printParserResults(parser)
sss = parser.ODEcalcString(scale = scale)
print 'ODEcalcString:'
print '------------------------------------------------'
print sss
print '------------------------------------------------'
cc = compile(sss, 'bof.log','exec')
exec cc
print 'Result from exec:', calcDerivs
    
print '\n================================================'
print ' Example from TSH2a.txt'
print '================================================\n'

print "Parameters are"
for z in zip(par_names, m_Parameters):
    print "%-8s= %f" % z
print
print 'At time =',t[0]
for z in zip(var_names, y0):
    print "%-8s= %f" % z

print '\n------------------------------------------------'

reps = 200

from scipy import integrate
print 'Integrating...(repeating %d times)' % reps

begin = time.time()


# Stiff problem. Requires analytic Jacobian.

r = integrate.ode(calcDerivs).set_integrator('vode', method='bdf')
msg = 'OK'
Y = empty((data.shape[0], data.shape[1]-1))

itmax = len(times)
for i in xrange(reps):
    y0 = copy(data[0, 1:])
    Y[0] = y0
    r.set_initial_value(y0)
    r.t = times[0]
    it = 1
    while it < itmax:
        y1 = r.integrate(times[it])
        Y[it] = y1
        if not r.successful():
            msg = "Something wrong"
            break
        it +=1


end = time.time()
print 'Elapsed time: %.3f s' % (end-begin)


print '------------------------------------------------'
print '\nDone:',
print msg
print "%d time points computed." % Y.shape[0]
print
print 'At time =',t[-1]
for z in zip(var_names, Y[-1]):
    print "%-8s= %f" % z

if msg == 'OK':
    print '------------------------------------------------'
    print "Score (computed %d times):" %reps
    begin = time.time()
    for i in xrange(reps):
        S = (Y- data[:, 1:])**2
        score = nansum(S)
    print score
    end = time.time()
    print 'Elapsed time: %.3f s' % (end-begin)
    print '------------------------------------------------'

#plot results...

SDLTSH, HTA = Y.T
f1 = p.figure()
p.plot(t, SDLTSH, 'r-', label='SDLTSH')
p.plot(t, HTA  , 'b-', label='HTA')
p.grid()
p.legend(loc='best')
p.xlabel('time (s)')
p.ylabel('concentrations (mM)')
p.title('Example from time course TSH2a.txt')
p.show()
