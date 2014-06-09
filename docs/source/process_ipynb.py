import sys
import subprocess
import shutil
import os

nbfile = 'solving.ipynb'
exe_nbfile = 'solving_executed.ipynb'

print '--------------------------------------------'
print 'executing notebook', nbfile
aaa = ['runipy', nbfile, exe_nbfile]
print subprocess.check_output(aaa)

print '--------------------------------------------'
print 'converting notebook', nbfile
aaa = 'ipython nbconvert --to rst --template non_exec_rst.tpl solving.ipynb'.split()
print subprocess.check_output(aaa)
print 'converting notebook', exe_nbfile
aaa = 'ipython nbconvert --to rst --template non_exec_rst.tpl solving_executed.ipynb'.split()
print subprocess.check_output(aaa)

print '--------------------------------------------'
exe_nbfile = exe_nbfile.replace('.ipynb', '.rst')
print 'fixing img links...'
print 'reading', exe_nbfile

with open (exe_nbfile) as f:
    alltext = f.read()

alltext = alltext.replace('.. image:: solving_executed_files%5C','.. image:: _static/solving/')

print 'writing', exe_nbfile
with open (exe_nbfile, 'w') as f:
    f.write(alltext)


print '--------------------------------------------'
print 'copying files to _static dir...'
#copy solving_executed_files\*.* .\_static\solving
files = os.listdir( './solving_executed_files' )
for f in files:
    f = os.path.join('./solving_executed_files',f)
    print 'copying', f
    shutil.copy(f, '.\_static\solving')
print '--------- Done -----------------------------'