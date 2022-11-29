"""
Functions relating to particle generation
"""
import numpy
from warnings import warn
from at.physics import ohmi_envelope, radiation_parameters
from at.constants import clight
from at.lattice import AtError, AtWarning

__all__ = ['beam', 'sigma_matrix', 'sig_matrix']


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


def _sigma_matrix_lattice(ring=None, twiss_in=None, emitx=None, emity=None,
                          blength=None, espread=None, verbose=False):

    if ring is None and (emitx is None or emity is None or blength is None
                         or espread is None or twiss_in is None):
        raise AtError('Provide either a ring or twiss_in and ALL emittance '
                      'paramaters: (ex, ey, espread, blength)')

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
    sig_mat = _sigma_matrix_from_R66(rmat6, emitx, emity, blength, espread)
    return sig_mat


def sig_matrix(ring=None, **kwargs):
    """
    Calculate the correlation matrix to be used for particle generation

    Keyword Parameters:
        ring:           Lattice object or list of twiss parameters.
        twiss_in:       Data structure containing input twiss parameters.
        betax:          Input horizontal beta function [m]
        alphax:         Input horizontal alpha function [m]
        emitx:          Horizontal emittance [m.rad]
        betay:          Input vertical beta function [m]
        alphay:         Input vertical alpha function [m]
        emity:          Vertical emittance [m.rad]
        blength:        One sigma bunch length [m]
        espread:        One sigma energy spread [dp/p]
        verbose:        Boolean flag on whether to print information
                        to the terminal
    Returns:
        sigma_matrix:    6x6 correlation matrix

    If *twiss_in* is provided, its 'R' or 'alpha', 'beta' fields are used,
    horizontal and vertical emittances must be provided.

    Otherwise, if *ring* is provided, horizontal and vertical emittances must
    be provided if ring is 4d. For a 6d ring, the provided emittances overload
    the computed ones.

    Otherwise, *alphax*, *betax*, *emitxx*, *alphay*, *betaxy* and *emitxx*
    must be provided

    If the lattice object is provided ``ohmi_envelope`` is used to
    compute the correlated sigma matrix and missing emittances
    and longitudinal parameters

    If emittances and longitudinal parameters are provided, then
    the 2x2 uncorrelated matrices are computed for each plane (x,y,z)
    using the initial optics computed from ring.get_optics,
    and are combined into the 6x6 matrix.

    If neither a lattice object nor a twiss_in is provided,
    then the ``beta``, ``alpha`` and ``emittance`` for horizontal and
    vertical are required, as well as ``blength`` and ``espread``.
    This then computes the analytical uncoupled sigma matrix
    """
    def require(message, *keys):
        miss = [key for key in keys if key not in kwargs]
        if len(miss) == 0:
            return [kwargs[key] for key in keys]
        else:
            raise AtError('{0}: missing {1}'.format(message, ', '.join(miss)))

    def ignore(message, *keys):
        ok = [key for key in keys if key in kwargs]
        if any(ok):
            warn(AtWarning('{0}: {1} ignored'.format(message, ', '.join(ok))))

    def r_4d(ax, ay, bx, by):
        r6 = numpy.zeros((3, 6, 6))
        r6[0, :2, :2] = numpy.array([[bx, -ax], [-ax, (1+ax**2)/bx]])
        r6[1, 2:4, 2:4] = numpy.array([[by, -ay], [-ay, (1+ay**2)/by]])
        return r6

    def expand(rmat):
        if rmat.shape[0] < 3:
            r6 = numpy.zeros((3, 6, 6))
            r6[:2, :4, :4] = rmat
        else:
            r6 = rmat
        return r6

    if ring is not None:
        kwargs['ring'] = ring
    emit = None
    if 'twiss_in' in kwargs:
        ignore('"twiss_in" is provided', 'ring', 'alphax', 'alphay',
               'betax', 'betay')
        twin = kwargs['twiss_in']
        try:
            Rin = twin['R']
        except (ValueError, KeyError):  # record arrays throw ValueError !
            Rm = r_4d(twin['alpha'][0], twin['beta'][0],
                      twin['alpha'][1], twin['beta'][1])
        else:
            Rm = expand(Rin)
    elif 'ring' in kwargs:
        ignore('"ring" is provided', 'alphax', 'alphay', 'betax', 'betay')
        ring = kwargs['ring']
        if ring.is_6d:
            _, beam0, _ = ohmi_envelope(ring)
            Rm = beam0.mode_matrices
            emit = beam0.mode_emittances
        else:
            el0, _, _ = ring.linopt6()
            Rm = expand(el0['R'])
    else:
        alphax, alphay, betax, betay = require(
            'Neither "twiss_in" nor "ring" is provided',
            'alphax', 'alphay', 'betax', 'betay')
        Rm = r_4d(alphax, alphay, betax, betay)

    if emit is None:
        emx, emy = require('Emittances must be provided', 'emitx', 'emity')
        ems = None
    else:
        emx = kwargs.pop('emitx', emit[0])
        emy = kwargs.pop('emity', emit[1])
        ems = emit[2]
        if any(('espread' in kwargs, 'blength' in kwargs)):
            blength, espread = require('Both "blength" and "espread" required',
                                       'blength', 'espread')
            ems = blength * espread

    sigm = emx * Rm[0, :, :] + emy * Rm[1, :, :]
    if ems is None:
        if any(('espread' in kwargs, 'blength' in kwargs)):
            blength, espread = require('Both "blength" and "espread" required',
                                       'blength', 'espread')
            Rm[2, 4:6, 4:6] = numpy.array([[espread/blength, 0],
                                           [0, blength/espread]])
            ems = blength*espread
            sigm += ems * Rm[2, :, :]
        else:
            warn(AtWarning('4D beam'))

    return sigm


def sigma_matrix(ring=None, twiss_in=None, betax=None, alphax=None,
                 emitx=None, betay=None, alphay=None, emity=None,
                 blength=None, espread=None, verbose=False):
    """
    Calculate the correlation matrix to be used for particle generation

    Parameters:
        ring:           Lattice object or list of twiss parameters.
        twiss_in:       Data structure containing input twiss parameters.
        betax:          Input horizontal beta function [m]
        alphax:         Input horizontal alpha function [m]
        emitx:          Horizontal emittance [m.rad]
        betay:          Input vertical beta function [m]
        alphay:         Input vertical alpha function [m]
        emity:          Vertical emittance [m.rad]
        blength:        One sigma bunch length [m]
        espread:        One sigma energy spread [dp/p]
        verbose:        Boolean flag on whether to print information
                        to the terminal
    Returns:
        sigma_matrix:    6x6 correlation matrix

    If the lattice object is provided ``ohmi_envelope`` is used to
    compute the correlated sigma matrix and missing emittances
    and longitudinal parameters

    If emittances and longitudinal parameters are provided, then
    the 2x2 uncorrelated matrices are computed for each plane (x,y,z)
    using the initial optics computed from ring.get_optics,
    and are combined into the 6x6 matrix.

    If ``twiss_in`` is provided, it has to be produced with ``linopt6`` and
    contain the rmatrix. ``ring`` is used only to compute missing emittances.

    If neither a lattice object nor a twiss_in is provided,
    then the ``beta``, ``alpha`` and ``emittance`` for horizontal and
    vertical are required, as well as ``blength`` and ``espread``.
    This then computes the analytical uncoupled sigma matrix
    """
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


def beam(nparts: int, sigma, orbit=None):
    """
    Generates an array of random particles according to the given sigma
    matrix

    Parameters:
        nparts:         Number of particles
        sigma:          sigma_matrix as calculated by at.sigma_matrix
        orbit:          An orbit can be provided to give a center of
                        mass offset to the distribution
    Returns:
        particle_dist:  a (6, ``nparts``) matrix of coordinates
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
