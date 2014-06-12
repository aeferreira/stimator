import sys
import subprocess
import shutil
import os

nbfiles = 'solving.ipynb'

def process_list(nbfiles):
    exe_nbfiles = nbfiles.replace('.ipynb', '_executed.ipynb')
    print ('--------- executing notebook {0} ------------'.format(nbfiles))
    aaa = ['runipy', nbfiles, exe_nbfiles]
    print (subprocess.check_output(aaa))

    for n in [nbfiles, exe_nbfiles]:
        print ('--------- converting notebook {0} -----------'.format(n))
        aaa = ('ipython nbconvert --to rst --template non_exec_rst.tpl %s'%n).split()
        print (subprocess.check_output(aaa))

    rst_exe_nbfiles = exe_nbfiles.replace('.ipynb', '.rst')
    print ('--------- fixing img links in {0} ---------'.format(rst_exe_nbfiles))

    print ('reading %s'% rst_exe_nbfiles)
    with open (rst_exe_nbfiles) as f:
        alltext = f.read()

    alltext = alltext.replace('.. image:: solving_executed_files%5C','.. image:: _static/solving/')

    print ('writing %s'% rst_exe_nbfiles)
    with open (rst_exe_nbfiles, 'w') as f:
        f.write(alltext)

    print ('------------- copying files to static dir -----------')
    ddest = '_static\solving'
    try:
        os.makedirs(ddest)
    except Exception:
        pass
    files = os.listdir( './solving_executed_files' )
    for f in files:
        f = os.path.join('./solving_executed_files',f)
        print 'copying', f, 'to', ddest
        shutil.copy(f, ddest)
    print ('--------- Done -----------------------------')

if __name__ == '__main__':
    process_list(nbfiles)