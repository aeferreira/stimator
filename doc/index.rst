.. S-timator documentation master file, created by
   sphinx-quickstart on Tue Jan 20 11:35:02 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

S-timator: dynamical systems modelling in python
================================================

``S-timator`` is a Python library to analyse **ODE-based models**
(also known as *dynamic* or *kinetic* models). These models are often found
in many scientific fields, particularly in Physics, Chemistry, Biology and
Engineering.

Some of the features that ``S-timator`` offers are

- **A mini language used to describe models**: models can be input as plain text 
  following a very simple and human-readable language.
- **Basic analysis**: numerical solution of ODE's, parameter scanning.
- **Parameter estimation** and **model selection**: given experimental data in
  the form of time series and constrains on model operating ranges,
  built-in numerical optimizers can find parameter values and assist you in the
  experimental design for model selection.

``S-timator`` is in an **alpha stage**: many new features will be available soon.

For a brief introduction to the package, you can read the :ref:`Basic use <introduction>`.

Much more detail can be found in the :ref:`tutorial <tutorial>`.

You can also browse the :ref:`API reference <api_ref>` to see the kind of tools that are available.

   
Contents:
---------

.. toctree::
   :maxdepth: 1

   introduction
   installing
   api

Tutorial:
---------

.. toctree::
   :maxdepth: 2

   basic_features
   models
   solving
   par_estimation

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

