import sys

fixfilename = sys.argv[1]
print 'reading', fixfilename

with open (fixfilename) as f:
    alltext = f.read()

alltext = alltext.replace('.. image:: solving_executed_files%5C','.. image:: _static/solving/')

print 'writing', fixfilename
with open (fixfilename, 'w') as f:
    f.write(alltext)

