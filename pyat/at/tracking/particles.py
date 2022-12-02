"""
Functions relating to particle generation
"""
import numpy as np
from warnings import warn
from ..physics import Orbit, ohmi_envelope
from ..lattice import AtError, AtWarning, Lattice, random

__all__ = ['beam', 'sigma_matrix']


# noinspection PyPep8Naming
def sigma_matrix(ring: Lattice = None, **kwargs):
    r"""
    Calculate the correlation matrix :math:`\Sigma` to be used for particle
    generation

    Parameters:
        ring:           Lattice object or list of twiss parameters.

    Keyword Args:
        beam:               Particle coordinates
        twiss_in:           Data structure containing input twiss parameters.
        betax (float):      Input horizontal beta function [m]
        alphax (float):     Input horizontal alpha function [m]
        emitx (float):      Horizontal emittance [m.rad]
        betay (float):      Input vertical beta function [m]
        alphay (float):     Input vertical alpha function [m]
        emity (float):      Vertical emittance [m.rad]
        blength (float):    One sigma bunch length [m]
        espread (float):    One sigma energy spread [dp/p]

    Returns:
        sigma_matrix:       (6, 6) :math:`\Sigma`-matrix

    *beam* is a (6, nparticles) coordinates array of a particle population.
    If provided, all other parameters are ignored.

    If *twiss_in* is provided, its *R* or *alpha* and *beta* fields are used,
    *emitx* and *emity* must be provided.

    If *ring* is provided:

    * for a 6d lattice, :py:func:`.ohmi_envelope` is called to compute the
      equilibrium emittances. *emitx*, *emity* may overload these values.
    * for a 4d lattice, *emitx*, *emity* must be provided.

    Otherwise, *alphax*, *betax*, *emitx*, *alphay*, *betay* and *emity*
    must be provided. The result is an uncoupled :math:`\Sigma`-matrix

    If *espread* and *blength* are not specified, the results is a
    monochromatic beam.
    """
    def require(mess, *keys):
        miss = [key for key in keys if key not in kwargs]
        if len(miss) == 0:
            return [kwargs.pop(key) for key in keys]
        else:
            raise AtError('{0}: missing {1}'.format(mess, ', '.join(miss)))

    def ignore(mess):
        ok = kwargs.keys()
        if any(ok):
            warn(AtWarning('{0}: {1} ignored'.format(mess, ', '.join(ok))))

    def r_4d(ax, ay, bx, by):
        r6 = np.zeros((3, 6, 6))
        r6[0, :2, :2] = np.array([[bx, -ax], [-ax, (1+ax**2)/bx]])
        r6[1, 2:4, 2:4] = np.array([[by, -ay], [-ay, (1+ay**2)/by]])
        return r6

    def expand(rmat):
        if rmat.shape[0] < 3:
            r6 = np.zeros((3, 6, 6))
            r6[:2, :4, :4] = rmat
        else:
            r6 = rmat
        return r6

    def sigma_from_rmat(rmat, emit):

        if emit is None:
            emx, emy = require('Emittances must be provided', 'emitx', 'emity')
            ems = None
        else:
            emx = kwargs.pop('emitx', emit[0])
            emy = kwargs.pop('emity', emit[1])
            ems = emit[2]
            if any(('espread' in kwargs, 'blength' in kwargs)):
                blength, espread = require(
                    'Both "blength" and "espread" are required',
                    'blength', 'espread')
                ems = blength * espread

        sigm = emx * rmat[0, :, :] + emy * rmat[1, :, :]
        if ems is None:
            if any(('espread' in kwargs, 'blength' in kwargs)):
                blength, espread = require(
                    'Both "blength" and "espread" are required',
                    'blength', 'espread')
                rmat[2, 4:6, 4:6] = np.array([[espread / blength, 0],
                                              [0, blength / espread]])
                ems = blength * espread
                sigm += ems * rmat[2, :, :]
            else:
                warn(AtWarning('Monochromatic beam: no energy spread'))
        else:
            sigm += ems * rmat[2, :, :]
        return sigm

    if ring is not None:
        kwargs['ring'] = ring

    if 'beam' in kwargs:
        message = "'beam' is provided"
        beam = kwargs.pop('beam')
        orbit = kwargs.pop('orbit', None)
        if orbit is None:
            orbit = np.mean(beam, axis=1)
        beam -= orbit.reshape((6, 1))
        sigmat = (beam @ beam.T) / beam.shape[1]

    elif 'twiss_in' in kwargs:
        twin = kwargs.pop('twiss_in')
        message = "'twiss_in' is provided"
#       orbit = kwargs.pop('orbit', twin.getattr('closed_orbit', np.zeros(6)))
        try:
            Rin = twin['R']
        except (ValueError, KeyError):  # record arrays throw ValueError !
            Rm = r_4d(twin['alpha'][0], twin['beta'][0],
                      twin['alpha'][1], twin['beta'][1])
        else:
            Rm = expand(Rin)
        sigmat = sigma_from_rmat(Rm, None)

    elif 'ring' in kwargs:
        ring = kwargs.pop('ring')
        message = "'ring' is provided"
        el0, _, _ = ring.linopt6(orbit=kwargs.pop('orbit', None))
        orbit = el0.closed_orbit
        Rm = expand(el0['R'])
        if ring.is_6d:
            _, beam0, _ = ohmi_envelope(ring, orbit=orbit)
            sigmat = sigma_from_rmat(Rm, beam0.mode_emittances)
        else:
            sigmat = sigma_from_rmat(Rm, None)

    else:
        message = "Uncoupled beam"
        alphax, alphay, betax, betay = require(
            'Neither "twiss_in" nor "ring" is provided',
            'alphax', 'alphay', 'betax', 'betay')
        # orbit = kwargs.pop('orbit', np.zeros(6))
        Rm = r_4d(alphax, alphay, betax, betay)
        sigmat = sigma_from_rmat(Rm, None)

    ignore(message)
    return sigmat


def beam(nparts: int, sigma, orbit: Orbit = None):
    r"""
    Generates an array of random particles according to the given
    :math:`\Sigma`-matrix

    Parameters:
        nparts:         Number of particles
        sigma:          :math:`\Sigma`-matrix as calculated by
          :py:func:`sigma_matrix`
        orbit:          An orbit can be provided to give a center of
          mass offset to the distribution

    Returns:
        particle_dist:  a (6, *nparts*) matrix of coordinates
    """
    def _get_single_plane(slc):
        try:
            return np.linalg.cholesky(sigma[slc, slc])
        except np.linalg.LinAlgError:
            return np.zeros((2, 2))

    v = random.thread.normal(size=(sigma.shape[1], nparts))

    try:
        lmat = np.linalg.cholesky(sigma)
    except np.linalg.LinAlgError:
        lmat = np.zeros((6, 6))
        lmat[:2, :2] = _get_single_plane(slice(2))
        lmat[2:4, 2:4] = _get_single_plane(slice(2, 4))
        lmat[4:, 4:] = _get_single_plane(slice(4, 6))

    particle_dist = np.asfortranarray(lmat @ v)

    if orbit is not None:
        particle_dist += orbit.reshape((6, 1))

    return particle_dist
