Welcome to AT's documentation!
================================


Accelerator Toolbox is a code used for simulating particle accelerators, used
particularly for synchrotron light sources. It is hosted on `Github
<https://github.com/atcollab/at>`_. Its original implementation is in Matlab.

pyAT is a Python interface to Accelerator Toolbox. It uses the 'pass methods'
defined in Accelerator Toolbox, implemented by compiling the C code used in the
AT 'integrators' into a Python extension. These pass methods are used by
higher-level functions to provide physics results.

.. toctree::
   :maxdepth: 2

   Installation
   primer
   demo1

.. toctree::
   :maxdepth: 2
   :caption: How to:
   :hidden:

   howtos/CavityControl

.. toctree::
   :maxdepth: 2
   :caption: Release notes:
   :hidden:

   releases/r2.3
