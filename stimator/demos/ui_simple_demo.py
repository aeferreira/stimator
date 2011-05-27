# first try
print dir()
import time
t0 = time.time()
print 'O----------'

for i in range(10):
    print i

## print 'zzzzzzzzzzzzzzzzzz for 5.0 s'
## time.sleep(5.0)

ui.checkpoint()

print 
for i in range(50):
    for j in [1,2,3]:
        print 'Hi',i+1, j

ui.checkpoint()
print '!----------'

tf = time.time()
print 'done in %f s'%(tf-t0)

res = ui.ok_cancel('Are you OK?', 'Status')
if res:
    print "you're OK"
else:
    print "you're NOT OK"
ui.demo_plot()
figure = ui.new_figure()
figure.add_subplot(1,1,1).plot([0,1,2], [2,4,6], 'ro')
