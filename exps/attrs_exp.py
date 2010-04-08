
class MyModel(object):
    pass

class MyDict(object):
    pass

m = MyModel()

m.V = 2.0
m.init = {'x1':3.0, 'x2':1.5}
m.init2 = MyDict()
for x in m.init:
    setattr(m.init2, x, m.init[x])

print m.V
print m.init
print m.init2

setattr(m, 'V2' , 1.2)

print m.V2

setattr(m.init2, 'x3' , 7.7)

print m.init2.x3
