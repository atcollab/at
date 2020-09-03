"""
Functions relating to particle generation
"""
import numpy
from at.physics import ohmi_envelope
from at.lattice.lattice_object import Lattice

__all__ = ['beam', 'sigma_matrix']


def sigma_matrix(argin):
    """
    Calculate the correlation matrix to be used for particle generation

    PARAMETERS
        argin           Lattice object or list of twiss parameters.
                        Depending on the length of the list the function
                        will do different things:
                        [espread, blength] 
                            computes 2x2 longitudinal sigma matrix
                        [beta, alpha, emit] 
                            computes 2x2 transverse sigma matrix
                        [betax, alphax, emitx, betay, alphay, emity] 
                            computes the 4x4 sigma matrix for horizontal 
                            and vertical.
                        [betax, alphax, emitx, betay, alphay, emity, espread, 
                            blength] 
                            computes the 6x6 sigma matrix.

    KEYWORDS
        twiss=False     Flag which states whether the input is a
                        lattice object or list of twiss parameters

    OUTPUT
        sigma_matrix    correlation matrix (either 2x2, 4x4 or 6x6)

    """
    if isinstance(argin, Lattice):
        argin = argin.radiation_on(copy=True)

        emit0, beamdata, emit = ohmi_envelope(argin, refpts=[0])
        sig_matrix = emit.r66[0]


    else:
        if len(argin) == 2:
            [espread, blength] = argin
            sig_matrix = numpy.array([
                            [espread*espread, 0],
                            [0, blength*blength]
                                     ])

        elif len(argin) == 3:
            [bx, ax, epsx] = argin
            sig_matrix = epsx * numpy.array([
                                    [bx, -ax],
                                    [-ax, (1 + ax * ax)/bx]
                                            ])

        elif len(argin) == 6:
            [bx, ax, epsx, by, ay, epsy] = argin
            sig_matrix = numpy.block([
                            [sigma_matrix([bx, ax, epsx]),
                                numpy.zeros((2, 2))],
                            [numpy.zeros((2, 2)),
                                sigma_matrix([by, ay, epsy])]
                                    ])

        elif len(argin) == 8:
            [bx, ax, epsx, by, ay, epsy, espread, blength] = argin
            sig_matrix = numpy.block([
                            [sigma_matrix([bx, ax, epsx]),
                                numpy.zeros((2, 4))],
                            [numpy.zeros((2, 2)),
                                sigma_matrix([by, ay, epsy]),
                                numpy.zeros((2, 2))],
                            [numpy.zeros((2, 4)),
                                sigma_matrix([espread, blength])]
                                    ])
        else:
            raise AttributeError('Wrong number of inputs provided')

    return sig_matrix


def beam(np, sigma, orbit=None):
    """
    Generates an array of random particles according to the given sigma
    matrix

    PARAMETERS
        np              Number of particles
        sigma           sigma_matrix as calculated by at.sigma_matrix

    KEYWORDS
        orbit=None      An orbit can be provided to give a center of 
                        mass offset to the distribution

    OUTPUT
        particle_dist   a matrix of shape (M, np) where M is shape of
                        sigma matrix
    """
    v = numpy.random.normal(size=(sigma.shape[0], np))

    try:
        lmat = numpy.linalg.cholesky(sigma)

    except numpy.linalg.LinAlgError:
        row_idx = numpy.array([0, 1, 4, 5])
        a1 = numpy.linalg.cholesky(sigma[row_idx[:, None], row_idx])
        a = numpy.block([[a1, numpy.zeros((4, 2))], [numpy.zeros((2, 6))]])
        row_idx = numpy.array([0, 1, 4, 5, 2, 3])
        lmat = a[row_idx[:, None], row_idx]

    particle_dist = numpy.squeeze(numpy.dot(lmat, v))

    if orbit is not None:
        if (not isinstance(orbit, (numpy.ndarray, list)) or
                len(orbit) != sigma.shape[0]):
            raise AttributeError('orbit should be a list or array' +
                                 ' with a length the same as sigma' +
                                 ' matrix')
        else:
            for i, orb in enumerate(orbit):
                particle_dist[i, :] += orb

    return particle_dist
