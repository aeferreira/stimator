import sys
import subprocess
import shutil
import os

nbfiles = [{'nb':'solving.ipynb', 'name':'solving'}]

def process_list(nbfiles):
    for nbf in nbfiles:
        name = nbf['nb']
        
        print ('--------- backing up notebook {0} ------------'.format(name))
        shutil.copyfile(os.path.realpath(name), os.path.realpath(name)+".bak")
        
        print ('--------- executing notebook {0} ------------'.format(name))
        aaa = ['runipy', '-o', name]
        print (subprocess.check_output(aaa))

        print ('--------- converting notebook {0} -----------'.format(name))
        aaa = ('ipython nbconvert --to rst %s'% name).split()
        print (subprocess.check_output(aaa))

        # MS-Windows fix:
        rst_name = name.replace('.ipynb', '.rst')
        print ('--------- fixing img links in {0} ---------'.format(rst_name))

        print ('reading %s'% rst_name)
        with open (rst_name) as f:
            alltext = f.read()

        alltext = alltext.replace('_files%5C', '_files/')

        print ('writing %s' % rst_name)
        with open (rst_name, 'w') as f:
            f.write(alltext)

        print ('--------- restoring non-executed {0} -------'.format(name))
        os.remove(name)
        os.rename(os.path.realpath(name)+".bak", os.path.realpath(name))

        print ('--------- Done -----------------------------')

if __name__ == '__main__':
    process_list(nbfiles)