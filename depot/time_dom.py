import numpy
import time

def dominance(vec1, vec2):
    """Compute Pareto dominance relationship."""
    d_result = 0
    for vo,vn in zip(vec1, vec2):
        d = vn-vo
        if d <= 0 and d_result <=0:
            d_result -= 1
        elif d >= 0 and d_result >=0:
            d_result += 1
        else:
            d_result = 0
            break
    return d_result

def numdominance(vec1, vec2):
    vec1 = numpy.array(vec1)
    vec2 = numpy.array(vec2)
    size = len(vec1)
    d = vec2 <= vec1
    if numpy.all(d): return -size
    d = vec2 >= vec1
    if numpy.all(d): return size
    return 0

## v1 = numpy.array([-1.0,-2.0, -3, -1, -2, -3])
## v2 = numpy.array([0.0,-1.0, -3, 0, -2, 0])
## v3 = numpy.array([0.0,-2.0, -3, -4.0, -2, -3])
v1 = [-1.0,-2.0, -3, -1, -2, -3]
v2 = [0.0,-1.0, -3, 0, -2, 0]
v3 = [0.0,-2.0, -3, -4.0, -2, -3]

print '---------------------dominance() ----------------------'
print 'dominance(v1, v2) = ', dominance (v1, v2)
print 'dominance(v2, v1) = ', dominance (v2, v1)
print 'dominance(v1, v3) = ', dominance (v1, v3)
print 'dominance(v3, v1) = ', dominance (v3, v1)

print '---------------------numnumdominance() ----------------------'
print 'numdominance(v1, v2) = ', numdominance (v1, v2)
print 'numdominance(v2, v1) = ', numdominance (v2, v1)
print 'numdominance(v1, v3) = ', numdominance (v1, v3)
print 'numdominance(v3, v1) = ', numdominance (v3, v1)
print '------------------- timings -------------------------------------'
repeats = 1000000
t0 = time.time()
dlist = []
for i in range(repeats):
    dlist.append(dominance(v2,v1))
elapsed = time.time() -t0
print 'dominance() took', elapsed

t0 = time.time()
ndlist = []
for i in range(repeats):
    ndlist.append(numdominance(v2,v1))
elapsed = time.time() -t0
print 'numdominance() took', elapsed

