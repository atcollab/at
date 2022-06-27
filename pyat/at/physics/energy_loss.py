from enum import Enum
from warnings import warn
from math import pi
from typing import Optional, Tuple
import numpy
from scipy.optimize import least_squares
from at.lattice import Lattice, Dipole, Wiggler, RFCavity, Refpts
from at.lattice import check_radiation, AtError
from at.lattice import checktype, set_value_refpts, get_cells, refpts_len
from at.constants import clight, Cgamma
from at.tracking import lattice_pass

__all__ = ['get_energy_loss', 'set_cavity_phase', 'ELossMethod',
           'get_timelag_fromU0']


class ELossMethod(Enum):
    """methods for the computation of energy losses"""
    #: The losses are obtained from
    #: :math:`E_{loss}=C_\gamma/2\pi . E^4 . I_2`.
    #: Takes into account bending magnets and wigglers.
    INTEGRAL = 1
    #: The losses are obtained by tracking without cavities.
    #: Needs radiation ON, takes into account all radiating elements
    TRACKING = 2


def get_energy_loss(ring: Lattice,
                    method: Optional[ELossMethod] = ELossMethod.INTEGRAL
                    ) -> float:
    """Computes the energy loss per turn

    Parameters:
        ring:           Lattice description
        method:         Method for energy loss computation.
          See :py:class:`ELossMethod`.

    Returns:
        eloss (float):  Energy loss per turn [eV]
    """

    # noinspection PyShadowingNames
    def integral(ring):
        """Losses = Cgamma / 2pi * EGeV^4 * i2
        """

        def wiggler_i2(wiggler: Wiggler):
            rhoinv = wiggler.Bmax / ring.BRho
            coefh = wiggler.By[1, :]
            coefv = wiggler.Bx[1, :]
            return wiggler.Length * (numpy.sum(coefh * coefh) + numpy.sum(
                coefv*coefv)) * rhoinv ** 2 / 2

        def dipole_i2(dipole: Dipole):
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
def get_timelag_fromU0(ring: Lattice,
                       method: ELossMethod = ELossMethod.INTEGRAL,
                       cavpts: Optional[Refpts] = None) -> Tuple[float, float]:
    """
    Get the TimeLag attribute of RF cavities based on frequency,
    voltage and energy loss per turn, so that the synchronous phase is zero.
    An error occurs if all cavities do not have the same frequency.
    Used in set_cavity_phase()

    Parameters:
        ring:               Lattice description
        method:             Method for energy loss computation.
          See :py:class:`ELossMethod`.
        cavpts:             Cavity location. If None, use all cavities.
          This allows to ignore harmonic cavities.
    Returns:
        timelag (float):    Timelag
        ts (float):         Time difference with the present value
    """
    def singlev(values):
        vals = numpy.unique(values)
        if len(vals) > 1:
            raise AtError('values not equal for all cavities')
        return vals[0]

    def eq(x):
        return numpy.sum(rfv*numpy.sin(2*pi*freq*(x+tl0)/clight))-u0

    def deq(x):
        return numpy.sum(rfv*numpy.cos(2*pi*freq*(x+tl0)/clight))

    if cavpts is None:
        cavpts = get_cells(ring, checktype(RFCavity))
    u0 = get_energy_loss(ring, method=method) / ring.periodicity
    freq = numpy.array([cav.Frequency for cav in ring.select(cavpts)])
    rfv = numpy.array([cav.Voltage for cav in ring.select(cavpts)])
    tl0 = numpy.array([cav.TimeLag for cav in ring.select(cavpts)])
    try:
        frf = singlev(freq)
        tml = singlev(tl0)
    except AtError:
        if u0 > numpy.sum(rfv):
            raise AtError('Not enough RF voltage: unstable ring')
        ctmax = 1/numpy.amin(freq)*clight/2
        tt0 = tl0[numpy.argmin(freq)]
        zero_diff = least_squares(deq, -tt0+ctmax/2,
                                  bounds=(-tt0, ctmax-tt0)).x[0]
        if numpy.sign(deq(zero_diff-1.0e-6)) > 0:
            ts = least_squares(eq, (zero_diff-tt0)/2,
                               bounds=(-tt0, zero_diff)).x[0]
        else:
            ts = least_squares(eq, (ctmax-tt0+zero_diff)/2,
                               bounds=(zero_diff, ctmax-tt0)).x[0]
        timelag = ts+tl0
    else:
        vrf = numpy.sum(rfv)
        if u0 > vrf:
            raise AtError('Not enough RF voltage: unstable ring')
        timelag = clight/(2*pi*frf)*numpy.arcsin(u0/vrf)
        ts = timelag - tml
        timelag *= numpy.ones(refpts_len(ring, cavpts))
    return timelag, ts


def set_cavity_phase(ring: Lattice,
                     method: ELossMethod = ELossMethod.TRACKING,
                     refpts: Optional[Refpts] = None,
                     cavpts: Optional[Refpts] = None,
                     copy: bool = False) -> None:
    """
    Adjust the TimeLag attribute of RF cavities based on frequency,
    voltage and energy loss per turn, so that the synchronous phase is zero.
    An error occurs if all cavities do not have the same frequency.

    .. Warning::

       This function changes the time reference, this should be avoided

    Parameters:
        ring:       Lattice description
        method:     Method for energy loss computation.
          See :py:class:`ELossMethod`.
        cavpts:     Cavity location. If None, use all cavities.
          This allows to ignore harmonic cavities.
        copy:       If True, returns a shallow copy of ring with new
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
