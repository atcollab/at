"""
Functions relating to particle generation
"""
import numpy
from at.physics import ohmi_envelope, radiation_parameters
from at.lattice.constants import clight
from at.lattice import AtError

__all__ = ['beam', 'sigma_matrix']


def _generate_2d_trans_matrix(emit, beta, alpha):
    return emit*numpy.array([[beta, -alpha], [-alpha, (1+alpha**2)/beta]])


def _generate_2d_long_matrix(espread, blength):
    return numpy.array([[espread*espread, 0], [0, blength*blength]])


def _generate_2d_long_Rmatrix(espread, blength):
    return numpy.array([[espread/blength, 0], [0, blength/espread]])


def _sigma_matrix_uncoupled(betax, alphax, emitx, betay, alphay, emity,
                            blength, espread):
    sig_matrix = numpy.zeros((6, 6))
    sig_matrix[:2, :2] = _generate_2d_trans_matrix(emitx, betax, alphax)
    sig_matrix[2:4, 2:4] = _generate_2d_trans_matrix(emity, betay, alphay)
    sig_matrix[4:, 4:] = _generate_2d_long_matrix(espread, blength)
    return sig_matrix


def _sigma_matrix_from_R66(R66, emitx, emity, blength, espread):
    sig_matrix = emitx*R66[0] + emity*R66[1] + \
             blength*espread*R66[2]
    return sig_matrix


def _compute_bunch_length_from_espread(ring, espread):
    rp = radiation_parameters(ring.radiation_off(copy=True))
    blength = ring.beta*clight*numpy.abs(rp.etac)*espread / \
        (2.0*numpy.pi*rp.f_s)
    return blength


def _sigma_matrix_lattice(ring, twiss_in=None, emitx=None, emity=None,
                          blength=None, espread=None, verbose=False):
    if espread is not None and blength is None:
        blength = _compute_bunch_length_from_espread(ring, espread)
    flag_all = emitx and emity and espread
    flag_any = emitx or emity or espread

    if not flag_all:
        if verbose:
            print('Calculating missing parameters (ex, ey or espread) '
                  'using ohmi envelop: needs radiations ON')
        try:
            emit0, beamdata, emit = ohmi_envelope(ring, refpts=[0])
        except AtError:
            raise AtError('Please provide ex, ey, espread or turn on '
                          'radiations to compute the sigma matrix')
        if emitx is None:
            emitx = beamdata.mode_emittances[0]
        if emity is None:
            emity = beamdata.mode_emittances[1]
        if espread is None:
            espread = numpy.sqrt(emit0.r66[4, 4])
            blength = numpy.sqrt(emit0.r66[5, 5])

    if not flag_any and not twiss_in:
        return emit.r66[0]
    elif twiss_in:
        if verbose:
            print('Generating pseudo-correlated matrix '
                  'from twiss_in')
        if not hasattr(twiss_in, 'R'):
            raise AtError('twiss_in has to the contain the R matrix. '
                          'Please use the output from linopt6.')
        rmat = twiss_in.R
    else:
        if verbose:
            print('Generating pseudo-correlated matrix '
                  'from start point of ring')
        l0, _, _ = ring.get_optics()
        rmat = l0.R

    if rmat.shape[0] != 3:
        rmat6 = numpy.zeros((3, 6, 6))
        rmat6[0, :4, :4] = rmat[0]
        rmat6[1, :4, :4] = rmat[1]
        rmat6[2, 4:, 4:] = _generate_2d_long_Rmatrix(espread, blength)
    else:
        rmat6 = rmat
    sig_matrix = _sigma_matrix_from_R66(rmat6, emitx, emity, blength,
                                        espread)
    return sig_matrix


def sigma_matrix(ring=None, twiss_in=None, betax=None, alphax=None,
                 emitx=None, betay=None, alphay=None, emity=None,
                 blength=None, espread=None, verbose=False):
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
    if ring is not None:
        return _sigma_matrix_lattice(ring, twiss_in=twiss_in,
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
        lmat = numpy.zeros((6, 6))
        lmat[:2, :2] = _get_single_plane([0, 1])
        lmat[2:4, 2:4] = _get_single_plane([2, 3])
        lmat[4:, 4:] = _get_single_plane([4, 5])

    particle_dist = numpy.squeeze(numpy.dot(lmat, v))
    if particle_dist.ndim == 1:
        particle_dist = numpy.array([particle_dist]).T

    if orbit is not None:
        if numpy.shape(orbit) != (6,):
            raise AtError('beam: input orbit shape has to be (6,)')
        particle_dist = (particle_dist.T + numpy.array(orbit)).T

    return particle_dist
