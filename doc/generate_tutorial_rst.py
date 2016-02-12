import sys
import subprocess
import shutil
import os
import glob

nbfiles = [{'nb':'basic_features.ipynb', 'name':'basic_features'},
           {'nb':'models.ipynb', 'name':'models'},
           {'nb':'solving.ipynb', 'name':'solving'},
           {'nb':'par_estimation.ipynb', 'name':'par_estimation'}]

def process_list(nbfiles):
    print ('=========== Assembling supporting files ======')
    imgs = glob.glob('../notebooks/images/*.*')
    for img in imgs:
        print ('-- copying file {0}'.format(img))
        shutil.copy(img, './images')
    
    for nbf in nbfiles:
        name = nbf['nb']
        
        print ('=========== notebook {0} ================'.format(name))
        print ('-- copying')
        fromname = "../notebooks/" + name
        shutil.copy(fromname, '.')
        
        print ('-- executing')
        aaa = ('ipython nbconvert --execute --inplace --to notebook %s'% name).split()
        print (subprocess.check_output(aaa))

        rst_name = name.replace('.ipynb', '.rst')
        print ('-- converting to {0}'.format(rst_name))
        aaa = ('ipython nbconvert --to rst %s'% name).split()
        print (subprocess.check_output(aaa))

        # MS-Windows fix:
        print ('-- fixing img links in {0}'.format(rst_name))

        print ('reading %s'% rst_name)
        with open (rst_name) as f:
            alltext = f.read()

        alltext = alltext.replace('_files%5C', '_files/')
        alltext = alltext.replace("    %matplotlib inline\n", '')

        print ('writing %s' % rst_name)
        with open (rst_name, 'w') as f:
            f.write(alltext)

    print ('=========== Done =====================')

if __name__ == '__main__':
    process_list(nbfiles)