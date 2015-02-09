# -*- coding: utf8 -*-
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import re
import os

here = os.path.abspath(os.path.dirname(__file__))

def find_version(*file_paths):
    p = os.path.join(here, *file_paths)
    with open(p) as f:
        version_file = f.read()
    version_match = re.search(r"^__version__\.version = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

with open('README.rst') as f:
    readme = f.read()
## with open('HISTORY.rst') as f:
##     history = f.read()

packages = [
    'stimator', 
    'stimator.gui',
    'stimator.tests',
    'stimator.examples',
    'stimator.moo',
    'stimator.smallseaborn'
]

requires = ['sympy', 'pandas', 'seaborn']

classifs=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Natural Language :: English',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 2 :: Only',
        'Topic :: Scientific/Engineering :: Artificial Life',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics']
    
setup(name = "stimator",
    version=find_version('stimator', '__init__.py'),
    license = "BSD",
    description = "Analysis of ODE models with focus on model selection and parameter estimation.",
    author = "António Ferreira",
    author_email = "aeferreira@fc.ul.pt",
    url = "http://enzymology.fc.ul.pt/software.htm",
    include_package_data=True,
    packages = packages,
    entry_points = {
        'gui_scripts': [
            'start_stimator_wx = stimator.gui:run_wxgui',
        ]
    },
    keywords = "ODE-models estimation dynamics",
    classifiers=classifs,
    long_description = readme,
    install_requires = requires)


