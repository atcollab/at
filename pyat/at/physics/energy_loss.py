from enum import Enum
from warnings import warn
from math import pi
import numpy
from scipy.optimize import least_squares
from at.lattice import Lattice, Dipole, Wiggler, RFCavity
from at.lattice import check_radiation, AtError, AtWarning
from at.lattice import checktype, set_value_refpts, get_cells, refpts_len
from at.lattice.constants import clight, Cgamma
from at.tracking import lattice_pass

__all__ = ['get_energy_loss', 'set_cavity_phase', 'ELossMethod',
           'get_timelag_fromU0']


class ELossMethod(Enum):
    """Enum class for energy loss methods"""
    INTEGRAL = 1
    TRACKING = 2


def get_energy_loss(ring, method=ELossMethod.INTEGRAL):
    """Compute the energy loss per turn [eV]

    PARAMETERS
        ring                        lattice description

    KEYWORDS
        method=ELossMethod.INTEGRAL method for energy loss computation
            The enum class ELossMethod declares 2 values
            INTEGRAL: The losses are obtained from
                Losses = Cgamma / 2pi * EGeV^4 * i2
                Takes into account bending magnets and wigglers.
            TRACKING: The losses are obtained by tracking without cavities.
                Needs radiation ON, takes into account all radiating elements.
    """

    # noinspection PyShadowingNames
    def integral(ring):
        """Losses = Cgamma / 2pi * EGeV^4 * i2
        """

        def wiggler_i2(wiggler):
            rhoinv = wiggler.Bmax / ring.BRho
            coefh = wiggler.By[1, :]
            coefv = wiggler.Bx[1, :]
            return wiggler.Length * (numpy.sum(coefh * coefh) + numpy.sum(
                coefv*coefv)) * rhoinv ** 2 / 2

        def dipole_i2(dipole):
            return dipole.BendingAngle ** 2 / dipole.Length

        i2 = 0.0
        for el in ring:
            if isinstance(el, Dipole):
                i2 += dipole_i2(el)
            elif isinstance(el, Wiggler) and el.PassMethod != 'DriftPass':
                i2 += wiggler_i2(el)
        e_loss = Cgamma / 2.0 / pi * ring.energy ** 4 * i2
        return e_loss

    # noinspection PyShadowingNames
    @check_radiation(True)
    def tracking(ring):
        """Losses from tracking
        """
        ringtmp = ring.radiation_off(dipole_pass=None,
                                     quadrupole_pass=None,
                                     wiggler_pass=None,
                                     sextupole_pass=None,
                                     octupole_pass=None,
                                     copy=True)
        diff_elem = ringtmp.get_elements('Diffusion')
        for de in diff_elem:
            ringtmp.remove(de)

        o6 = numpy.squeeze(lattice_pass(ringtmp, numpy.zeros(6),
                           refpts=len(ringtmp)))
        if numpy.isnan(o6[0]):
            dp = 0
            for e in ringtmp:
                ot = numpy.squeeze(lattice_pass([e], numpy.zeros(6)))
                dp += -ot[4] * ring.energy
            return dp
        else:
            return -o6[4] * ring.energy

    if isinstance(method, str):
        method = ELossMethod[method.upper()]
        warn(FutureWarning('You should use {0!s}'.format(method)))
    if method is ELossMethod.INTEGRAL:
        return ring.periodicity * integral(ring)
    elif method == ELossMethod.TRACKING:
        return ring.periodicity * tracking(ring)
    else:
        raise AtError('Invalid method: {}'.format(method))


# noinspection PyPep8Naming
def get_timelag_fromU0(ring, method=ELossMethod.INTEGRAL, cavpts=None):
    """
    Get the TimeLag attribute of RF cavities based on frequency,
    voltage and energy loss per turn, so that the synchronous phase is zero.
    An error occurs if all cavities do not have the same frequency.
    Used in set_cavity_phase()

    PARAMETERS
        ring        lattice description

    KEYWORDS
        method=ELossMethod.TRACKING
                    method for energy loss computation. See "get_energy_loss".
        cavpts=None Cavity location. If None, use all cavities.
                    This allows to ignore harmonic cavities.
    RETURN
        timelag     Timelag
        ts          Time difference with the present value
    """
    def singlev(values):
        vals = numpy.unique(values)
        if len(vals) > 1:
            raise AtError('values not equal for all cavities')
        return vals[0]

    def eq(x, freq, rfv, tl0, u0):
        omf = 2*numpy.pi*freq/clight   
        eq1 = numpy.sum(-rfv*numpy.sin(omf*(x-tl0)))-u0
        eq2 = numpy.sum(-omf*rfv*numpy.cos(omf*(x-tl0)))
        if eq2 > 0:
            return numpy.sqrt(eq1**2+eq2**2)
        else:
            return eq1

    if cavpts is None:
        cavpts = ring.get_cells(checktype(RFCavity))
    u0 = ring.get_energy_loss(method=method) / ring.periodicity
    freq = numpy.array([cav.Frequency for cav in ring.select(cavpts)])
    rfv = numpy.array([cav.Voltage for cav in ring.select(cavpts)])
    tl0 = numpy.array([cav.TimeLag for cav in ring.select(cavpts)])
    if u0 > numpy.sum(rfv):
        raise AtError('Not enough RF voltage: unstable ring')
    try:
        frf = singlev(freq)
        tml = singlev(tl0)
    except AtError:
        ctmax = clight/numpy.amin(freq)/2
        tt0 = tl0[numpy.argmin(freq)]
        ts = ctmax/numpy.pi*numpy.arcsin(u0/numpy.sum(rfv))
        bounds = (-ctmax+tt0, ctmax+tt0)
        args = (freq,rfv,tl0,u0)
        r = [least_squares(eq,bounds[0]/4-ts, args=args, bounds=bounds),
             least_squares(eq,bounds[1]/4-ts, args=args, bounds=bounds)]
        res = numpy.array([abs(ri.fun[0]) for ri in r])
        ok = res < 1.0e-6
        if not numpy.any(ok):
            raise AtError('No solution found for Phis, please check '
                          'RF settings')
        if numpy.all(ok) and abs(r[0].x[0]-r[1].x[0])>1.0e-6:
            warn(AtWarning('More than one solution found for Phis: use '
                           'best fit, please check RF settings'))
        ts = -r[numpy.argmin(res)].x[0]
        timelag = ts+tl0
    else:
        vrf = numpy.sum(rfv)
        timelag = clight/(2*numpy.pi*frf)*numpy.arcsin(u0/vrf)
        ts = timelag - tml
        timelag *= numpy.ones(at.refpts_len(ring,cavpts))
    return timelag, ts


def set_cavity_phase(ring, method=ELossMethod.TRACKING,
                     refpts=None, cavpts=None, copy=False):
    """
   Adjust the TimeLag attribute of RF cavities based on frequency,
   voltage and energy loss per turn, so that the synchronous phase is zero.
   An error occurs if all cavities do not have the same frequency.

   !!!!WARNING!!!: This function changes the time reference,
   this should be avoided

    PARAMETERS
        ring        lattice description

    KEYWORDS
        method=ELossMethod.INTEGRAL
                            method for energy loss computation.
                            See "get_energy_loss".
        cavpts=None         Cavity location. If None, use all cavities.
                            This allows to ignore harmonic cavities.
        copy=False          If True, returns a shallow copy of ring with new
                            cavity elements. Otherwise, modify ring in-place.
    """
    # refpts is kept for backward compatibility
    if cavpts is None and refpts is not None:
        warn(FutureWarning('You should use "cavpts" instead of "refpts"'))
        cavpts = refpts
    elif cavpts is None:
        cavpts = get_cells(ring, checktype(RFCavity))
    timelag, _ = get_timelag_fromU0(ring, method=method, cavpts=cavpts)
    set_value_refpts(ring, cavpts, 'TimeLag', timelag, copy=copy)


Lattice.get_energy_loss = get_energy_loss
Lattice.energy_loss = property(get_energy_loss)
Lattice.set_cavity_phase = set_cavity_phase
