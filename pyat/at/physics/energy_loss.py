import sys
from enum import Enum
from warnings import warn
from math import sqrt, pi
import numpy
from at.lattice import Lattice, Dipole, Wiggler
from at.lattice import check_radiation, set_cavity, AtError
from at.tracking import lattice_pass
from at.physics import clight, Cgamma, e_mass

__all__ = ['get_energy_loss', 'set_cavity_phase', 'ELossMethod']


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
    def integral(ring):
        """Losses = Cgamma / 2pi * EGeV^4 * i2
        """

        def wiggler_i2(wiggler):
            rhoinv = wiggler.Bmax / Brho
            coefh = wiggler.By[1, :]
            coefv = wiggler.Bx[1, :]
            return wiggler.Length * (numpy.sum(coefh * coefh) + numpy.sum(
                coefv*coefv)) * rhoinv ** 2 / 2

        def dipole_i2(dipole):
            return dipole.BendingAngle ** 2 / dipole.Length

        Brho = sqrt(ring.energy**2 - e_mass**2) / clight
        i2 = 0.0
        for el in ring:
            if isinstance(el, Dipole):
                i2 += dipole_i2(el)
            elif isinstance(el, Wiggler) and el.PassMethod != 'DriftPass':
                i2 += wiggler_i2(el)
        e_loss = Cgamma / 2.0 / pi * ring.energy ** 4 * i2
        return e_loss

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
        o0 = numpy.zeros(6)
        o6 = numpy.squeeze(lattice_pass(ringtmp, o0, refpts=len(ringtmp)))
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


def set_cavity_phase(ring, method=ELossMethod.INTEGRAL,
                     refpts=None, cavpts=None, copy=False):
    """
   Adjust the TimeLag attribute of RF cavities based on frequency,
   voltage and energy loss per turn, so that the synchronous phase is zero.
   An error occurs if all cavities do not have the same frequency.

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
    rfv = ring.get_rf_voltage(cavpts=cavpts)
    freq = ring.get_rf_frequency(cavpts=cavpts)
    print("\nThis function modifies the time reference\n"
          "This should be avoided, you have been warned!\n",
          file=sys.stderr)
    u0 = get_energy_loss(ring, method=method)
    if u0 > rfv:
        raise AtError('Not enough RF voltage: unstable ring')
    timelag = clight / (2*pi*freq) * numpy.arcsin(u0/rfv)
    set_cavity(ring, TimeLag=timelag, cavpts=cavpts, copy=copy)


Lattice.get_energy_loss = get_energy_loss
Lattice.energy_loss = property(get_energy_loss)
Lattice.set_cavity_phase = set_cavity_phase
