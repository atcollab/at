from __future__ import print_function
import numpy
# noinspection PyUnresolvedReferences
from .atpass import atpass
from .lattice import uint32_refpts


def linepass(line, r_in, refpts=None, keep_lattice=False):
    """linepass tracks particles through each element of the sequence line
    calling the element-specific tracking function specified in the
    line[i].PassMethod field.

    r_out = linepass(line, r_in) tracks particle(s) with initial
    condition(s) r_in to the end of the line

    line        AT lattice
    r_in        6xN matrix: input coordinates of N particles

    r_out       6xN matrix: output coordinates of N particles at

    r_out=linepass(line, r_in, REFPTS) returns intermediate results
    at the entrance of each element specified in refpts

    REFPTS is an array of increasing indexes that selects elements
    between 0 and length(line).
    See further explanation of REFPTS in the 'help' for FINDSPOS
    r_out       6x(N*length(REFPTS)) matrix: output coordinates of N particles at
                each reference point

    keep_lattice:If set to True, linepass is more efficient because it reuses
                 some of the data and functions stored in the persistent
                 memory in previous calls to linepass.

    !!! In order to use this option, linepass must first be called without
    the keep_lattice flag. This will create persistent data structures and keep
    pointers to pass-method functions.

    !!! linepass(..., keep_lattice=True) assumes that the RING is unchanged
    since the last call. Otherwise, ringpass with keep_lattice=False must be
    called again.

    NOTE:
    linepass(line,r_in,length(line)) is the same as linepass(line,r_in)
    since the reference point length(line) is the exit of the last element
    linepass(line,r_in, 0) is a copy of r_in since the
    reference point 0 is the entrance of the first element
    """
    if refpts is None:
        refpts = len(line)
    refs = uint32_refpts(refpts, len(line))
    r_in = numpy.asfortranarray(r_in.reshape((6, -1)))
    return atpass(line, r_in, 1, refs, int(keep_lattice))


def ringpass(ring, r_in, nturns=1, keep_lattice=False):
    """ringpass tracks particles through each element of the cell array ring
    calling the element-specific tracking function specified in the
    ring[i].PassMethod field.

    r_out=ringpass(ring,r_in,nturns) tracks particle(s) with initial
    condition(s) RIN for nturns turns

    ring:       AT lattice
    r_in:       6xN matrix: input coordinates of N particles
    nturns:     Number of turns to perform (default: 1)

    r_out:      6x(N*nturns) matrix: output coordinates of N particles at

    keep_lattice:If set to True, ringpass is more efficient because it reuses
                 some of the data and functions stored in the persistent
                 memory in previous calls to ringpass.

    !!! In order to use this option, ringpass must first be called without
    the keep_lattice flag. This will create persistent data structures and keep
    pointers to pass-method functions.

    !!! ringpass(..., keep_lattice=True) assumes that the ring is unchanged
    since the last call. Otherwise, ringpass with keep_lattice=False must be
    called again.
    """
    r_in = numpy.asfortranarray(r_in.reshape((6, -1)))
    return atpass(ring, r_in, nturns, reuse=int(keep_lattice))
