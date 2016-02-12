import sys
import subprocess
import shutil
import os
import glob
from os import path as pth

nbfiles = [{'nb':'basic_features.ipynb', 'name':'basic_features'},
           {'nb':'models.ipynb', 'name':'models'},
           {'nb':'solving.ipynb', 'name':'solving'},
           {'nb':'par_estimation.ipynb', 'name':'par_estimation'}]

def process_list(nbfiles):
    imgs = glob.glob('images/*.*')
    for img in imgs:
        print ('-- deleting file {0}'.format(img))
        os.remove(img)
    
    for nbf in nbfiles:
        name = nbf['nb']
        
        if pth.isfile(name):
            print ('-- deleting notebook {0}'.format(name))
            os.remove(name)
        
        rst_name = name.replace('.ipynb', '.rst')
        if pth.isfile(rst_name):
            print ('-- deleting {0}'.format(rst_name))
            os.remove(rst_name)
        
        dir_name = name[:-6] + '_files'
        if pth.isdir(dir_name):
            print ('-- deleting dir {0}'.format(dir_name))
            shutil.rmtree(dir_name)

    print ('--------- Done')

if __name__ == '__main__':
    process_list(nbfiles)