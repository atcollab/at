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

   Installation <getting_started/01_Installation>
   AT Primer <getting_started/02_Primer>
   Example <getting_started/03_Lattice_example>

.. toctree::
   :maxdepth: 2
   :caption: How to:
   :hidden:

   Control the RF cavities <howtos/100_Control_RF_cavities>

.. toctree::
   :maxdepth: 2
   :caption: User guide
   :hidden:

   api/lattice
   api/atphysics
   api/atplot
   api/attrack
   api/atmatch
   api/atutils
   api/atmat

.. toctree::
   :maxdepth: 2
   :caption: Release notes:
   :hidden:

   release_notes/r2.3
   release_notes/r2.4
   release_notes/r2.5
   release_notes/r2.6
   release_notes/r2.7
