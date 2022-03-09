"""
Functions relating to particle generation
"""
import numpy
from at.physics import ohmi_envelope
from at.lattice.constants import clight, e_mass
from at.physics import get_mcf, get_tune
from at.lattice import AtError, AtWarning
import warnings

__all__ = ['beam', 'sigma_matrix']


def _generate_2d_trans_matrix(emit, beta, alpha):
    return emit*numpy.array([[beta, -alpha], [-alpha, (1+alpha**2)/beta]])


def _generate_2d_long_matrix(espread, blength):
    return numpy.array([[espread*espread, 0], [0, blength*blength]])


def _sigma_matrix_uncoupled(betax, alphax, emitx,
                            betay, alphay, emity,
                            blength, espread):

    sig_matrix_x = _generate_2d_trans_matrix(emitx,
                                             betax,
                                             alphax)
    sig_matrix_y = _generate_2d_trans_matrix(emity,
                                             betay,
                                             alphay)
    sig_matrix_long = _generate_2d_long_matrix(espread, blength)

    sig_matrix = numpy.block([
                             [sig_matrix_x, numpy.zeros((2, 4))],
                             [numpy.zeros((2, 2)), sig_matrix_y,
                              numpy.zeros((2, 2))],
                             [numpy.zeros((2, 4)), sig_matrix_long]
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


def _sigma_matrix_from_R66(R66, emitx, emity, blength, espread):
    sig_matrix = emitx*R66[0] + emity*R66[1] + \
             blength*espread*R66[2]
    return sig_matrix


def _sigma_matrix_lattice(ring=None, twiss_in=None, emitx=None,
                          emity=None, blength=None, espread=None,
                          verbose=False):

    flag = emitx or emity or espread
    if flag:
        assert emitx is not None, 'emitx must be defined'
        assert emity is not None, 'emity must be defined'
        assert espread is not None, 'espread must be defined'

    if twiss_in and not flag:
        raise AtError('Emittances must be specified for twiss_in')
    if twiss_in:
        if not hasattr(twiss_in, 'R'):
            raise AtError('twiss_in should contain the R matrix. '
                          'Please use the output from linopt6.')

    if flag and twiss_in:
        assert blength is not None, 'blength must be defined for twiss_in'
    elif flag and ring:
        if espread is not None and blength is None:
            blength = _compute_bunch_length_from_espread(ring, espread)

    if ring and not flag:

        cavPassFlag = numpy.any([i.PassMethod == 'CavityPass' for i in ring])
        radPassFlag = numpy.any(['Rad' in i.PassMethod for i in ring])

        if cavPassFlag and not radPassFlag and not flag:
            raise AtError('Cannot compute 6D sigma matrix without '
                          'radiation damping and without emittances. '
                          'Either switch on radiation or provide the '
                          'emittances.')

        if ring.radiation:
            if verbose:
                print('Generating correlated sigma matrix using '
                      'ohmi envelope')
            emit0, beamdata, emit = ohmi_envelope(ring, refpts=[0])
            sig_matrix = emit.r66[0]

    elif flag:
        if ring:
            ld0, bd, ld = ring.get_optics()
            rmat = ld0.R
        elif twiss_in:
            rmat = twiss_in.R

        if verbose:
            print('Generating pseudo-correlated matrix '
                  'from initial optics conditions')

        sig_matrix = _sigma_matrix_from_R66(rmat,
                                            emitx, emity, blength, espread)

    else:
        raise AtError('A lattice or twiss_in must be provided')
    return sig_matrix


def sigma_matrix(ring=None, twiss_in=None, **kwargs):

    """
    Calculate the correlation matrix to be used for particle generation

    PARAMETERS
        ring            Lattice object or list of
                        twiss parameters.
        twiss_in        Data structure containing input
                        twiss parameters.

    KEYWORDS
        betax           Input horizontal beta function [m]
        alphax          Input horizontal alpha function [m]
        emitx           Horizontal emittance [m.rad]
        betay           Input vertical beta function [m]
        alphay          Input vertical alpha function [m]
        emity           Vertical emittance [m.rad]
        blength         One sigma bunch length [m]
        espread         One sigma energy spread [dp/p]
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

    If neither a lattice object nor a twiss_in is provided,
    then the beta, alpha and emittance for horizontal and
    vertical is required, as well as blength and espread.
    This then computes the analytical uncoupled sigma matrix
    """

    d = kwargs
    betax = d.get('betax', None)
    alphax = d.get('betax', None)
    emitx = d.get('emitx', None)
    betay = d.get('betay', None)
    alphay = d.get('betay', None)
    emity = d.get('emity', None)
    blength = d.get('blength', None)
    espread = d.get('espread', None)
    verbose = d.get('verbose', False)

    if ring is not None and twiss_in is not None:
        warnings.warn(AtWarning('Both ring and twiss_in have been provided. '
                                'Ignoring twiss_in and taking the ring.'))

    if ring is not None or twiss_in is not None:
        return _sigma_matrix_lattice(ring=ring, twiss_in=twiss_in,
                                     emitx=emitx, emity=emity,
                                     blength=blength, espread=espread,
                                     verbose=verbose)
    else:
        ln = '{:s} must be defined for the uncoupled sigma_matrix'
        assert betax is not None, ln.format('betax')
        assert alphax is not None, ln.format('alphax')
        assert emitx is not None, ln.format('emitx')
        assert betay is not None, ln.format('betay')
        assert alphay is not None, ln.format('alphay')
        assert emity is not None, ln.format('emity')
        assert blength is not None, ln.format('blength')
        assert espread is not None, ln.format('espread')

        return _sigma_matrix_uncoupled(betax, alphax, emitx,
                                       betay, alphay, emity,
                                       blength, espread)


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
        particle_dist   a matrix of shape (6, np)
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
