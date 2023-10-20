pyAT Examples
=============

Please follow ``<at>/pyat/README.rst`` for installation instructions.
In these examples, we use the ring from `hmba.mat` found in ``<at>/pyat/machine_data``;
however, they will also work on the ``err.mat`` ring if you want to be more
representative of a real accelerator.

Initialisation:
---------------

- Start Python within the virtual environment where PyAT in installed and
  import ``at``::

    $ cd <at>
    $ source venv/bin/activate
    $ python
    Python 2.7.3 (default, Nov  9 2013, 21:59:00)
    [GCC 4.4.7 20120313 (Red Hat 4.4.7-3)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import at
    >>>

- Load a pyAT ring from a .mat file::

    >>> ring = at.load_lattice('pyat/machine_data/hmba.mat')

Basic Use:
----------

- The lattice object is a Python interator::

    >>> for element in ring:
    ...   print(element)

- Viewing the first element in the ring::

    >>> ring[0]
    RFCavity('RFC', 0.0, 187500.0, 352372212.4670127, 31, 6000000000.0, PassMethod='IdentityPass')

- Changing the pass method of an element::

    >>> ring[0].PassMethod = "CavityPass"

- Get the s position of the 10th element in the ring::

    >>> at.get_s_pos(ring, 10)
    array([3.4295565])

Calculating Physics Data:
-------------------------

- Return the linear optics data at the first 5 elements::

    >>> at.get_optics(ring, refpts=[0, 1, 2, 3, 4])

- Physics data functions can also be called as methods on the lattice::

    >>> ring.find_orbit(refpts=[0, 1, 2, 3, 4])

