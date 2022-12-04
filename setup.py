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
    version_match = re.search(r"__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


with open('README.rst') as f:
    readme = f.read()


packages = [
    'stimator',
    'stimator.moo',
    'stimator.examples',
]

requires = ['six', 'pytest', 'sympy', 'numpy', 'scipy', 'matplotlib>=2.0']

classifs = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Natural Language :: English',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Artificial Life',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics']

sdesc = "Analysis of ODE models with parameter estimation and model selection."

setup(name="stimator",
      version=find_version('stimator', '__init__.py'),
      license="BSD",
      description=sdesc,
      author="António Ferreira",
      author_email="aeferreira@fc.ul.pt",
      url="http://webpages.fc.ul.pt/~aeferreira/stimator",
      include_package_data=True,
      packages=packages,
      package_data={'stimator.examples':
                    ['timecourses/TSH2a.txt', 'timecourses/TSH2b.txt']},
      keywords="ODE-models estimation dynamics",
      classifiers=classifs,
      long_description=readme,
      install_requires=requires)
