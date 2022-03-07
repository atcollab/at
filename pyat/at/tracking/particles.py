"""
Functions relating to particle generation
"""
import numpy
from at.physics import ohmi_envelope
from at.lattice.constants import clight, e_mass
from at.physics import get_mcf, get_tune
from at.lattice import AtWarning, AtError
import warnings

__all__ = ['beam', 'sigma_matrix']


def _generate_2d_trans_matrix(emit, beta, alpha):
    return emit*numpy.array([[beta, -alpha], [-alpha, (1+alpha**2)/beta]])


def _generate_2d_long_matrix(espread, blength):
    return numpy.array([[espread*espread, 0], [0, blength*blength]])


def _generate_sigma_matrix(ld0, emitx, emity, blength, espread, radiation):
    if radiation:
        try:
            sig_matrix = emitx*ld0.R[0] + emity*ld0.R[1] + \
                         blength*espread*ld0.R[2]
        except AttributeError:
            warnings.warn(AtWarning('R matrices not found. '
                                    '(are you using linopt6?)'))
            warnings.warn(AtWarning('Creating uncoupled sigma matrix'))

            sig_matrix_x = _generate_2d_trans_matrix(emitx,
                                                     ld0.beta[0],
                                                     ld0.alpha[0])
            sig_matrix_y = _generate_2d_trans_matrix(emity,
                                                     ld0.beta[1],
                                                     ld0.alpha[1])
            sig_matrix_long = _generate_2d_long_matrix(espread, blength)
            sig_matrix = numpy.block([
                                     [sig_matrix_x, numpy.zeros((2, 4))],
                                     [numpy.zeros((2, 2)), sig_matrix_y,
                                      numpy.zeros((2, 2))],
                                     [numpy.zeros((2, 4)), sig_matrix_long]
                                     ])
    else:
        sig_matrix_long = _generate_2d_long_matrix(espread, blength)

        try:
            sig_matrix_trans = (emitx*ld0.R[0] + emity*ld0.R[1])[:4, :4]
        except AttributeError:
            warnings.warn(AtWarning('R matrices not found. '
                                    '(are you using linopt6?)'))
            warnings.warn(AtWarning('Creating uncoupled sigma matrix'))

            sig_matrix_x = _generate_2d_trans_matrix(emitx,
                                                     ld0.beta[0],
                                                     ld0.alpha[0])
            sig_matrix_y = _generate_2d_trans_matrix(emity,
                                                     ld0.beta[1],
                                                     ld0.alpha[1])
            sig_matrix_trans = numpy.block([
                                     [sig_matrix_x, numpy.zeros((2, 2))],
                                     [numpy.zeros((2, 2)), sig_matrix_y],
                                     ])

        sig_matrix = numpy.block([
                        [sig_matrix_trans,
                            numpy.zeros((4, 2))],
                        [numpy.zeros((2, 4)),
                            sig_matrix_long]
                                ])
    return sig_matrix


def _compute_bunch_length_from_espread(ring, espread):
    gamma = ring.energy/e_mass
    beta = numpy.sqrt(1.0-1.0/gamma/gamma)
    f0 = clight/ring.circumference

    mcf = get_mcf(ring.radiation_off(copy=True))
    if ring.radiation:
        qs3 = get_tune(ring)
    else:
        qs3 = get_tune(ring.radiation_on(copy=True))

    f_s = qs3[2]*f0

    blength = clight * beta * numpy.abs(mcf) * espread / (2 * numpy.pi * f_s)
    return blength


def sigma_matrix(ring=None, twiss_in=None, emitx=None,
                 emity=None, blength=None, espread=None):
    """
    Calculate the correlation matrix to be used for particle generation

    PARAMETERS
        ring            Lattice object or list of
                        twiss parameters.
        twiss_in        Data structure containing input
                        twiss parameters.
        emitx           Horizontal emittance [m.rad]
        emity           Vertical emittance [m.rad]
        blength         One sigma bunch length [m]
        espread         One sigma energy spread [dp/p]

    KEYWORDS
        verbose=False   Boolean flag on whether to print information
                        to the terminal
    OUTPUT
        sigma_matrix    6x6 correlation matrix

    If the lattice object is provided with no other
    arguments, ohmi_envelope is used to compute the
    correlated sigma matrix.

    If the lattice object and emittances and longitudinal
    parameters are provided, then the 2x2 uncorrelated
    matrices are computed for each plane (x,y,z) using
    the initial optics computed from ring.get_optics,
    and are combined together into the 6x6 matrix.

    If the twiss_in is provided alongside the emittances and
    longitudinal parameters, then the 2x2 uncorrelated
    matrices are computed for each plane and combined
    into the 6x6 matrix.
    """
    flag = emitx or emity or blength or espread
    if flag:
        assert emitx is not None, 'emitx must be defined'
        assert emity is not None, 'emity must be defined'
        assert blength is not None, 'blength must be defined'
        assert espread is not None, 'espread must be defined'

    if ring:
        cavPassFlag = numpy.any(
                        numpy.array(
                          [i.PassMethod == 'CavityPass' for i in ring]))
        radPassFlag = numpy.any(
                        numpy.array(
                          ['Rad' in i.PassMethod for i in ring]))

        if cavPassFlag and not radPassFlag and not flag:
            raise AtError('Cannot compute 6D sigma matrix without '
                          'radiation damping and without emittances')

        ld0, bd, ld = ring.get_optics()
        if ring.radiation and not flag:
            warnings.warn(AtWarning('Generating correlated sigma matrix using '
                                    'ohmi envelope'))

            emit0, beamdata, emit = ohmi_envelope(ring, refpts=[0])
            sig_matrix = emit.r66[0]
        elif flag:
            warnings.warn(AtWarning('Generating pseudo-correlated matrix '
                      'from initial optics conditions'))
            if ring.radiation:
                warnings.warn(AtWarning('Ignoring provided blength and '
                          'calculating it based on espread'))

                blength = _compute_bunch_length_from_espread(ring, espread)
            sig_matrix = _generate_sigma_matrix(ld0, emitx, emity,
                                                blength, espread,
                                                ring.radiation)
        else:
            raise AtError('Radiation is off but no '
                                 'emittances are specified')
    elif twiss_in:
        if not emitx:
            raise AtError('Emittances must be '
                                 'specified for twiss_in')

        warnings.warn(AtWarning('Generating un-correlated sigma matrix '
                  'from parameters in twiss_in'))
        sig_matrix = _generate_sigma_matrix(twiss_in, emitx,
                                            emity, blength,
                                            espread, False)

    else:
        raise AtError('A lattice or twiss_in must be provided')
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
    def _get_single_plane(dims):
        row_idx = numpy.array(dims)
        try:
            return numpy.linalg.cholesky(sigma[row_idx[:, None], row_idx])
        except numpy.linalg.LinAlgError:
            return numpy.zeros((2, 2))

    v = numpy.random.normal(size=(sigma.shape[0], nparts))

    try:
        lmat = numpy.linalg.cholesky(sigma)
    except numpy.linalg.LinAlgError:
        warnings.warn(AtWarning('Decomposition failed for 6x6 correlation matrix.'
              ' Computing 3 planes individually'))
        a1 = _get_single_plane([0, 1])
        a2 = _get_single_plane([2, 3])
        a3 = _get_single_plane([4, 5])
        lmat = numpy.block([[a1, numpy.zeros((2, 4))],
                            [numpy.zeros((2, 2)), a2,
                             numpy.zeros((2, 2))],
                            [numpy.zeros((2, 4)), a3]])

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
            for i, orb in enumerate(orbit):
                particle_dist[i, :] += orb

    return particle_dist
