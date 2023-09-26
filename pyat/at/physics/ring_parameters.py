from math import pi, sqrt, asin, cos
import numpy
from numpy import nan
from typing import Optional
from .radiation import get_radiation_integrals, ohmi_envelope
from .energy_loss import get_energy_loss, get_timelag_fromU0
from ..lattice import Lattice, Orbit
from ..constants import clight, Cgamma, Cq

__all__ = ['RingParameters', 'radiation_parameters', 'envelope_parameters']


class RingParameters(object):
    """Class for pretty printing the ring properties"""

    _props = {
        'tunes':            '              Frac. tunes: {0}',
        'tunes6':           '  Frac. tunes (6D motion): {0}',
        'fulltunes':        '                    Tunes: {0}',
        'chromaticities':   '           Chromaticities: {0}',
        'alphac':           ' Momentum compact. factor: {0:e}',
        'etac':             '              Slip factor: {0:e}',
        'E0':               '                   Energy: {0:e} eV',
        'U0':               '       Energy loss / turn: {0:e} eV',
        'i1':               ' Radiation integrals - I1: {0} m',
        'i2':               '                       I2: {0} m^-1',
        'i3':               '                       I3: {0} m^-2',
        'i4':               '                       I4: {0} m^-1',
        'i5':               '                       I5: {0} m^-1',
        'emittances':       '          Mode emittances: {0}',
        'J':                'Damping partition numbers: {0}',
        'Tau':              '            Damping times: {0} s',
        'sigma_e':          '            Energy spread: {0:g}',
        'sigma_l':          '             Bunch length: {0:g} m',
        'voltage':          '         Cavities voltage: {0} V',
        'phi_s':            '        Synchrotron phase: {0:g} rd',
        'f_s':              '    Synchrotron frequency: {0:g} Hz'
    }

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        vrs = vars(self).copy()
        vals = [(self._props[k], vrs.pop(k, None)) for k in self._props]
        # Predefined attributes
        lines = [f.format(v) for f, v in vals if v is not None]
        # Other attributes
        lines += ['{0:>25}: {1}'.format(k, getattr(self, k)) for k in vrs]
        return '\n'.join(lines)


# noinspection PyPep8Naming
def radiation_parameters(ring: Lattice, dp: Optional[float] = None,
                         params: Optional[RingParameters] = None,
                         **kwargs):
    r"""Compute ring parameters from the radiation integrals

    Valid for uncoupled lattices with no RF cavity or radiating element.

    Parameters:
        ring:       Lattice description.
        dp:         Momentum deviation.
        params:     :py:class:`.RingParameters` object to be updated.
          Default: create a new one

    Keyword Args:
        dct (float):        Path lengthening. If specified, ``dp`` is
          ignored and the off-momentum is deduced from the path lengthening.
        method (Callable):  Method for linear optics:

          :py:obj:`~.linear.linopt2`: no longitudinal motion, no H/V coupling,

          :py:obj:`~.linear.linopt6` (default): with or without longitudinal
          motion, normal mode analysis
        orbit (Orbit):      Avoids looking for the closed orbit if is
          already known ((6,) array)
        XYStep (float):     Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):     Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`

    Returns:
        params:             :py:class:`.RingParameters` object.

    **params** is a :py:class:`.RingParameters` object with the following
    attributes:

    ==================  ========================================
    **tunes**           (3,) fractional (H, V, Long.) tunes
    **fulltunes**       (3,) full tunes
    **chromaticities**  (2,) H, V Chromaticities
    **alphac**          Momentum compaction factor
    **etac**            Frequency slip factor
    **E0**              Energy [eV]
    **U0**              Energy loss / turn [eV]
    **i1**              Radiation integrals - :math:`I_1 \quad [m]`
    **i2**              :math:`I_2 \quad [m^{-1}]`
    **i3**              :math:`I_3 \quad [m^{-2}]`
    **i4**              :math:`I_4 \quad [m^{-1}]`
    **i5**              :math:`I_5 \quad [m^{-1}]`
    **emittances**      (3,) Mode emittances
    **J**               (3,) Damping partition numbers
    **Tau**             (3,) Damping times [s]
    **sigma_e**         Energy spread
    **sigma_l**         Bunch length [m]
    **voltage**         Total accelerating voltage [V]
    **phi_s**           Synchrotron phase [rad]
    **f_s**             Synchrotron frequency [Hz]
    ==================  ========================================
    """
    rp = RingParameters() if params is None else params
    _, ringdata, twiss = ring.get_optics(refpts=range(len(ring) + 1), dp=dp,
                                         get_chrom=True, **kwargs)
    rp.chromaticities = ringdata.chromaticity * ring.periodicity
    integs = get_radiation_integrals(ring, dp=dp, twiss=twiss)
    rp.i1, rp.i2, rp.i3, rp.i4, rp.i5 = numpy.array(integs) * ring.periodicity
    circumference = ring.circumference
    voltage = ring.rf_voltage
    E0 = ring.energy
    gamma = ring.gamma
    gamma2 = gamma * gamma
    beta = sqrt(1.0 - 1.0/gamma2)
    U0 = Cgamma / 2.0 / pi * E0**4 * rp.i2
    Jx = 1.0 - rp.i4/rp.i2
    Jz = 1.0
    Je = 2.0 + rp.i4/rp.i2
    damping_partition_numbers = numpy.array([Jx, Jz, Je])
    revolution_period = circumference / clight / beta
    ct = 2.0 * E0 / U0 * revolution_period
    rp.E0 = E0
    rp.U0 = U0
    emitx = Cq * gamma2 * rp.i5 / Jx / rp.i2
    alphac = rp.i1 / circumference
    etac = 1.0/gamma2 - alphac
    rp.phi_s = (pi - asin(U0 / voltage)) if U0 <= voltage else nan
    nus = sqrt(abs(etac * ring.harmonic_number *
               voltage * cos(rp.phi_s) / E0 / 2.0 / pi)) / beta
    rp.voltage = voltage
    rp.f_s = nus / revolution_period
    rp.Tau = ct / damping_partition_numbers
    rp.J = damping_partition_numbers
    rp.sigma_e = sqrt(Cq * gamma2 * rp.i3 / Je / rp.i2)
    rp.sigma_l = beta * abs(etac) * circumference / 2.0 / pi / nus * rp.sigma_e
    rp.emittances = numpy.array([emitx, nan, rp.sigma_e*rp.sigma_l])
    ringtunes, _ = numpy.modf(ring.periodicity * ringdata.tune)
    if len(ringtunes) < 3:
        rp.tunes = numpy.concatenate((ringtunes, (nus,)))
    else:
        rp.tunes6 = ringtunes
    rp.fulltunes = ring.periodicity * twiss[-1].mu / 2.0 / pi
    rp.alphac = alphac
    rp.etac = etac
    return rp


# noinspection PyPep8Naming
def envelope_parameters(ring: Lattice,
                        params: Optional[RingParameters] = None,
                        orbit: Orbit = None,
                        keep_lattice: bool = False) -> RingParameters:
    r"""Compute ring parameters from ohmi_envelope

    Parameters:
        ring:           Lattice description.
        params:         :py:class:`.RingParameters` object to be updated.
          Default: create a new one
        orbit:          Avoids looking for the closed orbit if it is
          already known ((6,) array)
        keep_lattice:   Assume no lattice change since the
          previous tracking.

    Returns:
        params:             :py:class:`.RingParameters` object.

    **params** is a :py:class:`.RingParameters` object with the following
    attributes:

    ==================  ========================================
    **tunes6**          (3,) fractional (H, V, Long.) tunes (6D motion)
    **emittances**      (3,) Mode emittances
    **J**               (3,) Damping partition numbers
    **Tau**             (3,) Damping times [s]
    **sigma_e**         Energy spread
    **sigma_l**         Bunch length [m]
    **voltage**         Total accelerating voltage [V]
    **phi_s**           Synchrotron phase [rad]
    **f_s**             Synchrotron frequency [Hz]
    ==================  ========================================
    """
    rp = RingParameters() if params is None else params
    emit0, beamdata, emit = ohmi_envelope(ring, orbit=orbit,
                                          keep_lattice=keep_lattice)
    voltage = ring.rf_voltage
    rp.E0 = ring.energy
    rp.U0 = get_energy_loss(ring)
    rev_freq = ring.revolution_frequency
    rp.Tau = 1.0 / rev_freq / beamdata.damping_rates / ring.periodicity
    alpha = 1.0 / rp.Tau
    rp.J = 4.0 * alpha / numpy.sum(alpha)
    rp.tunes6, _ = numpy.modf(ring.periodicity * beamdata.tunes)
    rp.phi_s = pi - numpy.arcsin(rp.U0 / voltage)
    rp.voltage = voltage
    rp.f_s = rp.tunes6[2] * rev_freq
    rp.emittances = beamdata.mode_emittances
    rp.sigma_e = sqrt(emit0.r66[4, 4])
    rp.sigma_l = sqrt(emit0.r66[5, 5])
    return rp


Lattice.radiation_parameters = radiation_parameters
Lattice.envelope_parameters = envelope_parameters
