# -*- coding: utf8 -*-
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.rst') as f:
    readme = f.read()
## with open('HISTORY.rst') as f:
##     history = f.read()

packages = [
    'stimator', 
    'stimator.gui',
    'stimator.tests',
    'stimator.examples',
    'stimator.moo'
]

requires = ['sympy']

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
    version = "0.9.85",
    license = "BSD",
    description = "Analysis of ODE models with focus on model selection and parameter estimation.",
    author = "António Ferreira",
    author_email = "aeferreira@fc.ul.pt",
    url = "http://enzymology.fc.ul.pt/software.htm",
    include_package_data=True,
    packages = packages,
    package_data={'stimator': ['examples/*.py', 
                               'examples/*.ipynb', 
                               'examples/*.mdl', 
                               'examples/*.txt']},
    scripts = ["run_wxgui.py"],
    keywords = "hello world example examples",
    classifiers=classifs,
    long_description = readme,
    install_requires = requires)


