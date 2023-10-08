.. pyAT documentation master file, created by
   sphinx-quickstart on Wed May 11 09:16:07 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AT
=============
Introduction
------------

Accelerator Toolbox (AT) is a collection of tools to model storage rings and beam transport lines.

With AT, it is possible to:

- **create and manipulate accelerator lattice elements**,
- **track particles through the lattice**, selecting the appropriate integrator to represent the physics
- **compute accelerator parameters and beam properties**, generating new scripts or taking advantage of the existing ones

AT is based on a 6-D modular tracking engine written in C/C++ for efficiency.
Lattice manipulation and computation of accelerator physics parameters are provided
by two interfaces:

- a :doc:`Matlab interface <m/index>`, available as a Matlab toolbox,
- a :doc:`python interface <p/index>`, available as a python package.

See :doc:`common/indx` for the global documentation.

.. toctree::
   :maxdepth: 2
   :hidden:

   AT basics <common/indx>
   Python <p/index>
   Matlab <m/index>
