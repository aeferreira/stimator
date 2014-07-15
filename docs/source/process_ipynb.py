import sys
import subprocess
import shutil
import os

nbfiles = ['solving.ipynb']

def process_list(nbfiles):
    for nbf in nbfiles:
        fname, ext = nbf.split('.')
        exe_nbf = nbf.replace('.ipynb', '_executed.ipynb')
        print ('--------- executing notebook {0} ------------'.format(nbf))
        aaa = ['runipy', nbf, exe_nbf]
        print (subprocess.check_output(aaa))

        for n in [nbf, exe_nbf]:
            print ('--------- converting notebook {0} -----------'.format(n))
            aaa = ('ipython nbconvert --to rst --template non_exec_rst.tpl %s'%n).split()
            print (subprocess.check_output(aaa))

        rst_exe_nbf = exe_nbf.replace('.ipynb', '.rst')
        print ('--------- fixing img links in {0} ---------'.format(rst_exe_nbf))

        print ('reading %s'% rst_exe_nbf)
        with open (rst_exe_nbf) as f:
            alltext = f.read()

        alltext = alltext.replace('.. image:: %s'% fname + '_executed_files%5C','.. image:: _static/%s/' % fname)

        print ('writing %s' % rst_exe_nbf)
        with open (rst_exe_nbf, 'w') as f:
            f.write(alltext)

        print ('------------- copying files to static dir -----------')
        ddest = '_static\%s' % fname
        try:
            os.makedirs(ddest)
        except Exception:
            pass
        files = os.listdir( './%s_executed_files' % fname )
        for f in files:
            f = os.path.join('./%s_executed_files' % fname,f)
            print 'copying', f, 'to', ddest
            shutil.copy(f, ddest)
        print ('--------- Done -----------------------------')

if __name__ == '__main__':
    process_list(nbfiles)