from numpy import *

a = array([1,2,3])
b = array([2,3,4])

print a
print b

d = dot(a,b)
print 'dot =', d

inn = inner(a,b)
print 'inner = ', inn

v = vdot(a,b)
print 'vdot = ', v

a = array(([1,2,3], [1,2,3]))
b = array(([2,3], [3,4], [4,5]))
print
print a
print b

d = dot(a,b)
print 'dot =', d

#~ inn = inner(a,b)
#~ print 'inner = ', inn

v = vdot(a,b)
print 'vdot = ', v


a = array(([1,2,3], [1,2,3]))
b = array(([1,2,3], [2,3,4]))
print
print a
print b

inn = inner(a,b)
print 'inner = ', inn



