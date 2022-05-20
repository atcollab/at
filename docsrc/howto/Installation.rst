Installation
============

pyAT is a Python interface to Accelerator Toolbox. It uses the ‘pass
methods’ defined in Accelerator Toolbox, implemented by compiling the C
code used in the AT ‘integrators’ into a Python extension. These pass
methods are used by higher-level functions to provide physics results.

Requirements
------------

pyAT supports Python 3.6 to 3.9.

Installation
------------

From PyPI
~~~~~~~~~

Install accelerator-toolbox from PyPI:

.. code:: shell

   $ pip install accelerator-toolbox

From GitHub
~~~~~~~~~~~

1. Install git on your computer.

2. Download the latest version of AT:
   ``shell $ git clone https://github.com/atcollab/at.git``
3. Go to the pyAT installation directory

   .. code:: shell

      $ cd at/pyat

4. Build and install

   .. code:: shell

      $ python setup.py install

Usage
-----

Example:

.. code:: python

   from at.load import load_mat
   from at.physics import linopt
   ring = load_mat('test_matlab/hmba.mat')
   linopt(ring, refpts=range(5))

For more examples of how to use pyAT, see ``pyat_examples.rst``.
