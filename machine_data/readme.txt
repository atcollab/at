This directory contains example lattice files using the creation
functions rather than the older style relying on global variables.
Each one can be called by writing e.g.
ring=FODO;


FODO.m                      FODO lattice from AT2.0 Primer document
soleil.m                    soleil lattice
thomx.m                     ThomX lattice
esrf.m                      esrf lattice
australian_synchrotron.m    Australian synchrotron lattice
dba.m                       dba lattice with creation functions

These scripts can also be used to create .mat files, from which the lattices
can be loaded in pyAT, e.g.

    >>> RING = esrf()
    >>> save('esrf.mat', 'RING')
