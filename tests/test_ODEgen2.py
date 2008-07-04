import sys
sys.path.append('..')
import modelparser
#import stimator_timecourse
from numpy import *
import pylab as p
import time

modelText = """
variables SDLTSH HTA

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


print 'rate strings:'
print '------------------------------------------------'
for k in parser.rates:
    print 'RATE',k['name'], ':'
    print 'raw:', k['rate']
    print 'new:', parser.rateCalcString(k['rate'])
print '------------------------------------------------'


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

ratebytecode = [compile(parser.rateCalcString(k['rate']), 'bof.log','eval') for k in parser.rates]

#~ variables = y0
#~ v = [eval(r) for r in ratebytecode]

#~ for i,k in enumerate(parser.rates):
    #~ print 'RATE',k['name'], ':', v[i]

N = zeros((len(parser.variables),len(parser.rates)), dtype=float)
for m, srow in enumerate(parser.stoichmatrixrows):
    for i,k in enumerate(parser.rates):
        if srow.has_key(k['name']):
            N[m,i] = scale*srow[k['name']]

#~ print 'N ='
#~ Nmat = mat(N)
#~ print Nmat

#~ v = mat(array(v))
#~ print 'v ='
#~ print v

#~ dXdt = v * (Nmat.T)

#~ print 'dXdt = N*v ='
#~ print dXdt

NT = N.transpose()
nvars = range(len(parser.variables))
v = empty(len(parser.variables))

def calcDerivs(variables, t):
    global v
    for i in nvars:
        v[i] = eval(ratebytecode[i])
    return dot(v,NT)

#~ print 'calcDerivs((0,0),0)'
#~ print calcDerivs((0,0),0)
#~ print 'calcDerivs(y0,0)'
#~ print calcDerivs(y0,0)

print '\n------------------------------------------------'

reps = 2000

from scipy import integrate
print 'Integrating...(repeating %d times)' % reps

begin = time.time()
for i in xrange(reps):
    y0 = copy(data[0, 1:])
    Y, infodict = integrate.odeint(calcDerivs, y0, times, full_output=True)
end = time.time()
print 'Elapsed time: %.3f s' % (end-begin)


print '------------------------------------------------'
print '\nDone:',
print infodict['message']
print "%d time points computed." % Y.shape[0]
print
print 'At time =',t[-1]
for z in zip(var_names, Y[-1]):
    print "%-8s= %f" % z
#print Y

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
