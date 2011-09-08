# -*- coding: latin1 -*-
from distutils.core import setup

setup(name = "stimator",
    version = "0.9.8",
    description = "S-timator : a Python package for the analysis of ODE models.",
    author = "Antonio Ferreira",
    author_email = "aeferreira@fc.ul.pt",
    url = "http://enzymology.fc.ul.pt/software.htm",
    packages = ['stimator', 'stimator.gui'],
    package_data={'stimator': ['stimator/examples/*.*']},
    scripts = ["stimator_gui.py"],
    long_description = """S-timator : a Python package for the analysis of ODE models.

Copyright 2005-2010 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows.""" 
) 

