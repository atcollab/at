pyAT
====

Introduction
------------

Accelerator Toolbox is a code used for simulating particle accelerators, used
particularly for synchrotron light sources. It is hosted on `Github
<https://github.com/atcollab/at>`_. Its original implementation is in Matlab.

pyAT is a Python interface to the Accelerator Toolbox. It uses the 'pass methods'
defined in Accelerator Toolbox, implemented by compiling the C code used in the
AT 'integrators' into a Python extension. These pass methods are used by
higher-level functions to provide physics results.

See the `pyAT website <https://atcollab.github.io/at/p/index.html>`_ for a
more detailed introduction.

pyAT supports Python 3.7 to 3.12.

Installation
------------

Install accelerator-toolbox from PyPI::

    $ pip install accelerator-toolbox

Usage
-----

Example usage::

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

Developer Notes
---------------

Developer notes are in ``developers.rst``.

