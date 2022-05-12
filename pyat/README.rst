pyAT
====

Introduction
------------

Accelerator Toolbox is a code used for simulating particle accelerators, used
particularly for synchrotron light sources. It is hosted on `Github
<https://github.com/atcollab/at>`_. Its original implementation is in Matlab.

pyAT is a Python interface to Accelerator Toolbox. It uses the 'pass methods'
defined in Accelerator Toolbox, implemented by compiling the C code used in the
AT 'integrators' into a Python extension. These pass methods are used by
higher-level functions to provide physics results.

See the `pyAT website <https://at.readthedocs.io/en/latest/>`_ for a
more detailed introduction.

pyAT supports Python 3.6 to 3.10.

Installation
------------

Install accelerator-toolbox from PyPI::

    pip install accelerator-toolbox

Usage
-----

Example usage::

    >>> from at.load import load_mat
    >>> from at.physics import get_optics
    >>> ring = load_mat('test_matlab/hmba.mat')
    >>> get_optics(ring, refpts=range(5))

For more examples of how to use pyAT, see ``pyat_examples.rst``.

Developer Notes
---------------

Developer notes are in ``developers.rst``.

