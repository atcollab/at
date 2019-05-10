pyAT Examples
=============

Please follow ``at/pyat/README.rst`` for installation instructions.
In these examples, we use the ring from `hmba.mat` found in ``test_matlab``;
however, they will also work on the ``err.mat`` ring if you want to be more
representative of a real accelerator.

Initialisation:
---------------

- Start Python inside the ``pyat`` directory::

    $ cd pyat
    $ source venv/bin/activate
    $ python
    Python 2.7.3 (default, Nov  9 2013, 21:59:00)
    [GCC 4.4.7 20120313 (Red Hat 4.4.7-3)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

- Load a pyAT ring from a .mat file::

    >>> import at
    >>> ring = at.load.load_mat('test_matlab/hmba.mat')

Basic Use:
----------

- Viewing the first element in the ring::

    >>> ring[0]
    RingParam('S28d', 6000000000.0, Periodicity=32)

- Changing the pass method of an element::

    >>> ring[1].PassMethod = "CavityPass"

- Get the s position of the 10th element in the ring::

    >>> at.lattice.get_s_pos(ring, 10)
    array([3.4295565])

- Creating a lattice object::

    >>> lattice = at.lattice.Lattice(ring)

Calculating Physics Data:
-------------------------

- Return the linear optics data at the first 5 elements::

    >>> at.physics.linopt(ring, refpts=[0, 1, 2, 3, 4])

- Physics data functions can also be called as methods on the lattice::

    >>> lattice.find_orbit4(refpts=[0, 1, 2, 3, 4])
