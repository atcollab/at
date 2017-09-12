from __future__ import print_function
import numpy
# noinspection PyUnresolvedReferences
from .atpass import atpass
from .lattice import uint32_refpts


# noinspection PyPep8Naming
def linepass(line, r_in, refpts=None, KeepLattice=False):
    """LINEPASS tracks particles through each element of the cell array LINE
    calling the element-specific tracking function specified in the
    LINE[i].PassMethod field.

    R_OUT=linepass(LINE, R_IN) tracks particle(s) with initial
    condition(s) r_in to the end of the line

    LINE        AT lattice
    R_IN        6xN matrix: input coordinates of N particles

    R_OUT       6xN matrix: output coordinates of N particles at

    R_OUT=linepass(LINE, R_IN, REFPTS) returns intermediate results
    at the entrance of each element specified in refpts

    REFPTS is an array of increasing indexes that selects elements
    between 0 and length(LINE).
    See further explanation of REFPTS in the 'help' for FINDSPOS
    R_OUT       6x(N*length(REFPTS)) matrix: output coordinates of N particles at
                each reference point

    KeepLattice:If set to True, LINEPASS is more efficient because it reuses
                some of the data and functions stored in the persistent
                memory in previous calls to LINEPASS.

    !!! In order to use this option, LINEPASS must first be called without
    the KeepLattice flag. This will create persistent data structures and keep
    pointers to pass-method functions.

    !!! LINEPASS(..., KeepLattice=True) assumes that the RING is unchanged
    since the last call. Otherwise, RINGPASS with KeepLattice=False must be
    called again.

    NOTE:
    linepass(LINE,R_IN,length(LINE)) is the same as LINEPASS(LINE,R_IN)
    since the reference point length(LINE) is the exit of the last element
    linepass(LINE,R_IN, 0) is a copy of R_IN since the
    reference point 0 is the entrance of the first element
    """
    if refpts is None:
        refpts = len(line)
    refs = uint32_refpts(refpts, len(line))
    r_in = numpy.asfortranarray(r_in.reshape((6, -1)))
    return atpass(line, r_in, 1, refs, int(KeepLattice))


# noinspection PyPep8Naming
def ringpass(ring, r_in, nturns=1, KeepLattice=False):
    """RINGPASS tracks particles through each element of the cell array RING
    calling the element-specific tracking function specified in the
    RING[i].PassMethod field.

    R_OUT=ringpass(RING,R_IN,NTURNS) tracks particle(s) with initial
    condition(s) RIN for NTURNS turns

    RING:       AT lattice
    R_IN:       6xN matrix: input coordinates of N particles
    NTURNS:     Number of turns to perform (default: 1)

    R_OUT:      6x(N*NTURNS) matrix: output coordinates of N particles at

    KeepLattice:If set to True, RINGPASS is more efficient because it reuses
                some of the data and functions stored in the persistent
                memory in previous calls to RINGPASS.

    !!! In order to use this option, RINGPASS must first be called without
    the KeepLattice flag. This will create persistent data structures and keep
    pointers to pass-method functions.

    !!! RINGPASS(..., KeepLattice=True) assumes that the RING is unchanged
    since the last call. Otherwise, RINGPASS with KeepLattice=False must be
    called again.
"""
    r_in = numpy.asfortranarray(r_in.reshape((6, -1)))
    return atpass(ring, r_in, nturns, reuse=int(KeepLattice))
