Installation
============

pyAT is a Python interface to Accelerator Toolbox. It uses the ‘pass
methods’ defined in Accelerator Toolbox, implemented by compiling the C
code used in the AT ‘integrators’ into a Python extension. These pass
methods are used by higher-level functions to provide physics results.

Requirements
------------

pyAT supports Python 3.7 to 3.12.

Installation
------------

From PyPI
~~~~~~~~~

PyPI maintains PyAT versions for Linux, MacOS and Windows, and python versions:

- CPython 3.7, 3.8, 3.9, 3.10, 3.11, 3.12
- PyPI does not support PyPy, because of the missing numpy and scipy libraries

*Installation with minimal dependencies:*
.........................................

::

   $ pip install accelerator-toolbox

This minimal installation disables all plotting functions. Note that the plotting
functions can be enabled later by simply installing `matplotlib <https://matplotlib.org>`_.

*Standard installation:*
........................

::

   $ pip install "accelerator-toolbox[plot]"

This installs in addition the `matplotlib <https://matplotlib.org>`_ package and its
dependencies, which enables all plotting functions.

*Optional dependencies:*
........................

In addition to ``[plot]``, other modifiers are available for specific purposes:

``"accelerator-toolbox[dev]"``
    Installs the test framework: `pytest <https://docs.pytest.org/en/stable/>`_,
    `pytest-cov <https://pypi.org/project/pytest-cov/>`_ and
    `flake8 <https://flake8.pycqa.org/en/latest/>`_. This allows to run locally
    the test sequence executed on GitHub.

``"accelerator-toolbox[doc]"``
    Installs `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ and its
    utilities, allowing the generate locally the HTML documentation.

These modifiers may be combined. Example: ``pip install "accelerator-toolbox[plot,dev]"``

From conda-forge
~~~~~~~~~~~~~~~~

conda-forge maintains PyAT versions for Linux, MacOS and Windows, and python versions:

- CPython 3.8, 3.9, 3.10, 3.11, 3.12
- PyPy 3.9

PyAT can be installed in a conda environment with::

   $ conda install -c conda-forge accelerator-toolbox

From source
~~~~~~~~~~~

1. Install git on your computer.

2. Download the latest version of AT::

    $ git clone https://github.com/atcollab/at.git

3. Go to the pyAT installation directory::

    $ cd at

4. Build and install::

    $ pip install .

Parallel computation
~~~~~~~~~~~~~~~~~~~~~
PyAT can be compiled for parallel processing. See :ref:`parallel`

Usage
-----

Example::

    >>> import at
    >>> ring = at.Lattice.load('machine_data/hmba.mat')
    >>> print(at.radiation_parameters(ring))
              Frac. tunes: [0.2099983  0.34001317 0.00349013]
                    Tunes: [76.2099983  27.34001317]
           Chromaticities: [5.73409894 3.91761206]
    Momentum compact. factor: 8.506669e-05
              Slip factor: -8.505944e-05
                   Energy: 6.000000e+09 eV
       Energy loss / turn: 2.526189e+06 eV
    Radiation integrals - I1: 0.07179435013387388 m
                       I2: 0.13844595446798158 m^-1
                       I3: 0.003357584058614851 m^-2
                       I4: -0.07375725030666251 m^-1
                       I5: 5.281495714523264e-07 m^-1
          Mode emittances: [1.3148797e-10           nan           nan]
    Damping partition numbers: [1.53275121 1.         1.46724879]
            Damping times: [0.00872477 0.0133729  0.00911427] s
            Energy spread: 0.000934463
             Bunch length: 0.0030591 m
         Cavities voltage: 6000000.0 V
        Synchrotron phase: 2.70701 rd
    Synchrotron frequency: 1239.74 Hz

For more examples of how to use pyAT, see ``pyat_examples.rst``.
