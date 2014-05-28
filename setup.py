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

    
setup(name = "stimator",
    version = "0.9.85",
    license = "BSD",
    description = "Analysis of ODE models with focus on model selection and parametrization.",
    author = "António Ferreira",
    author_email = "aeferreira@fc.ul.pt",
    url = "http://enzymology.fc.ul.pt/software.htm",
    include_package_data=True,
    packages = packages,
    package_data={'stimator': ['examples/*.*']},
    scripts = ["run_wxgui.py"],
    long_description = readme,
    install_requires = requires)


