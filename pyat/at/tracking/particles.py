"""
Functions relating to particle generation
"""
import numpy
from at.physics import ohmi_envelope
from at.lattice.lattice_object import Lattice

__all__ = ['beam', 'sigma_matrix']

def _generate_uncorrelated_sigma_matrix(bx, ax, emitx, by, ay, emity, blength, espread):
    sig_matrix_long = numpy.array([
                         [espread*espread, 0],
                         [0, blength*blength]
                                  ])

    sig_matrix_x = emitx * numpy.array([
                              [bx, -ax],
                              [-ax, (1 + ax * ax)/bx]
                                      ])

    sig_matrix_y = emity * numpy.array([
                              [by, -ay],
                              [-ay, (1 + ay * ay)/by]
                                       ])

    sig_matrix = numpy.block([
                    [sig_matrix_x,
                        numpy.zeros((2, 4))],
                    [numpy.zeros((2, 2)),
                        sig_matrix_y,
                        numpy.zeros((2, 2))],
                    [numpy.zeros((2, 4)),
                        sig_matrix_long]
                            ])
    return sig_matrix

def sigma_matrix(ring=None, twiss_in=None, emitx=None, emity=None, blength=None, espread=None):
    """
    Calculate the correlation matrix to be used for particle generation

    PARAMETERS
        ring            Lattice object or list of twiss parameters.
        twiss_in        Data structure containing inpu twiss parameters.
        emitx           Horizontal emittance [m.rad]
        emity           Vertical emittance [m.rad]
        blength         One sigma bunch length [m]
        espread         One sigma energy spread [dp/p]

    OUTPUT
        sigma_matrix    6x6 correlation matrix

    If the lattice object is provided with no other arguments, ohmi_envelope is used 
    to compute the correlated sigma matrix.

    If the lattice object and emittances and longitudinal parameters are provided, then
    the 2x2 uncorrelated matrices are computed for each plane (x,y,z) using the initial
    optics computed from ring.get_optics, and are combined together into the 6x6 matrix.

    If the twiss_in is provided alongside the emittances and longitudinal parameters, then
    the 2x2 uncorrelated matrices are computed for each plane and combined into the 6x6
    matrix. 

    """
    flag = emitx or emity or blength or espread
    if flag:
        assert emitx is not None, 'emitx must be defined'
        assert emity is not None, 'emity must be defined'
        assert blength is not None, 'blength must be defined'
        assert espread is not None, 'espread must be defined'

    if ring:
        if isinstance(ring, Lattice):
            if flag:
                print ('Generating uncorrelated sigma matrix by combining 2x2 matrices computed from initial optics conditions')
                ld0, bd, ld = ring.get_optics()
                bx, by = ld0.beta
                ax, ay = ld0.alpha
                sig_matrix = _generate_uncorrelated_sigma_matrix(bx, ax, emitx, by, ay, emity, blength, espread)
            else:
                print ('Generating correlated sigma matrix using ohmi envelope')
                ring = ring.radiation_on(copy=True)
                emit0, beamdata, emit = ohmi_envelope(ring, refpts=[0])
                sig_matrix = emit.r66[0]
        else:
            raise AttributeError('Input must be lattice object')

    elif twiss_in:

        bx, by = twiss_in.beta
        ax, ay = twiss_in.alpha

        if (blength and not espread) or (espread and not blength):
            raise AttributeError('Both blength and espread have to be provided')

        print ('Generating uncorrelated sigma matrix by combining 2x2 matrices computed from parameters in twiss_in')
        sig_matrix = _generate_uncorrelated_sigma_matrix(bx, ax, emitx, by, ay, emity, blength, espread)

    else:
        raise AttributeError('A lattice or twiss_in must be provided')
    return sig_matrix


def beam(nparts, sigma, orbit=None):
    """
    Generates an array of random particles according to the given sigma
    matrix

    PARAMETERS
        nparts          Number of particles
        sigma           sigma_matrix as calculated by at.sigma_matrix

    KEYWORDS
        orbit=None      An orbit can be provided to give a center of 
                        mass offset to the distribution

    OUTPUT
        particle_dist   a matrix of shape (M, np) where M is shape of
                        sigma matrix
    """
    v = numpy.random.normal(size=(sigma.shape[0], nparts))

    try:
        lmat = numpy.linalg.cholesky(sigma)

    except numpy.linalg.LinAlgError:
        print('Decomposition failed for 6x6 correlation matrix. Computing 3 planes individually')
        try:
            row_idx = numpy.array([0, 1])
            a1 = numpy.linalg.cholesky(sigma[row_idx[:, None], row_idx])
        except numpy.linalg.LinAlgError:
            a1 = numpy.zeros((2,2))

        try:
            row_idx = numpy.array([2, 3])
            a2 = numpy.linalg.cholesky(sigma[row_idx[:, None], row_idx])
        except numpy.linalg.LinAlgError:
            a2 = numpy.zeros((2,2))

        try:
            row_idx = numpy.array([4, 5])
            a3 = numpy.linalg.cholesky(sigma[row_idx[:, None], row_idx])
        except numpy.linalg.LinAlgError:
            a3 = numpy.zeros((2,2))

        lmat = numpy.block([[a1, numpy.zeros((2, 4))], [numpy.zeros((2, 2)), a2, numpy.zeros((2, 2))], [numpy.zeros((2,4)), a3]])

    particle_dist = numpy.squeeze(numpy.dot(lmat, v))
    if particle_dist.ndim == 1:
        particle_dist = numpy.array([particle_dist]).T

    if orbit is not None:
        if (not isinstance(orbit, (numpy.ndarray, list)) or
                len(orbit) != sigma.shape[0]):
            raise AttributeError('orbit should be a list or array' +
                                 ' with a length the same as sigma' +
                                 ' matrix')
        else:
            shp = particle_dist.shape
            for i, orb in enumerate(orbit):
                particle_dist[i, :] += orb

    return particle_dist
