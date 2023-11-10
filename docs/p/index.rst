.. pyAT documentation master file, created by
   sphinx-quickstart on Wed May 11 09:16:07 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyAT's documentation!
================================


Accelerator Toolbox is a code used for simulating particle accelerators, used
particularly for synchrotron light sources. It is hosted on `Github
<https://github.com/atcollab/at>`_. Its original implementation is in Matlab.

pyAT is a Python interface to Accelerator Toolbox. It uses the 'pass methods'
defined in Accelerator Toolbox, implemented by compiling the C code used in the
AT 'integrators' into a Python extension. These pass methods are used by
higher-level functions to provide physics results.

Sub-packages
------------

.. toctree::
   :maxdepth: 2
   :caption: Guides:
   :hidden:

   howto/Installation
   howto/Primer
   notebooks/test_parameters
   notebooks/test_observables
   examples

.. toctree::
   :maxdepth: 2
   :caption: How to:
   :hidden:

   howto/multiprocessing
   howto/CavityControl
   howto/Collective

.. autosummary::
   :toctree: api
   :caption: Packages:
   :recursive:

   at.lattice
   at.latticetools
   at.tracking
   at.physics
   at.load
   at.matching
   at.acceptance
   at.collective
   at.plot
   at.constants

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
