from math import pi, sqrt, asin, cos
import numpy
from numpy import nan
from scipy.constants import c
from at.lattice import Lattice
from at.physics import Cgamma, Cq, e_mass

__all__ = ['RingParameters', 'radiation_parameters', 'envelope_parameters']


class RingParameters(object):
    """Class for pretty printing the ring properties"""

    props = {
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
        'f_s':              '    Synchrotron frequency; {0:g} Hz'
    }

    def __init__(self, **kwargs):
        """Initialisation"""
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        vrs = vars(self).copy()
        vals = [(self.props[k], vrs.pop(k, None)) for k in self.props]
        # Predefined attributes
        lines = [f.format(v) for f, v in vals if v is not None]
        # Other attributes
        lines += ['{0:>25}: {1}'.format(k, getattr(self, k)) for k in vrs]
        return '\n'.join(lines)


# noinspection PyPep8Naming
def radiation_parameters(ring, dp=0.0, params=None):
    """Compute ring parameters from the radiation integrals. Valid for
    uncoupled lattices with no RF cavity or radiating element.

    INPUT
        ring            Lattice object.
        dp=0.0          momentum deviation.

    KEYWORD
        params=None     RingParam object to be updated.

    OUTPUT
        params          RingParam object. The computed attributes are,

            tunes           (3,) fractional (H, V, Long.) tunes
            fulltunes       (3,) full tunes
            chromaticities  (2,) H, V Chromaticities
            alphac          Momentum compaction factor
            etac            Frequency slip factor
            E0              Energy [eV]
            U0              nergy loss / turn [eV]
            i1              Radiation integrals - I1 [m]
            i2                                    I2 [m^-1]
            i3                                    I3 [m^-2]
            i4                                    I4 [m^-1]
            i5                                    I5 [m^-1]
            emittances      (3,) Mode emittances
            J               (3,) Damping partition numbers
            Tau             (3,) Damping times [s]
            sigma_e         Energy spread
            sigma_l         Bunch length [m]
            voltage         Total accelerating voltage [V]
            phi_s           Synchrotron phase [rad]
            f_s             Synchrotron frequency [Hz]
    """
    rp = RingParameters() if params is None else params
    _, tunes, chroms, twiss = ring.linopt(dp, range(len(ring) + 1),
                                          get_chrom=True, coupled=False)
    rp.chromaticities = chroms * ring.periodicity
    integs = ring.get_radiation_integrals(dp, twiss=twiss)
    rp.i1, rp.i2, rp.i3, rp.i4, rp.i5 = numpy.array(integs) * ring.periodicity
    circumference = ring.circumference
    revolution_period = circumference / c
    voltage = ring.voltage
    E0 = ring.energy
    gamma = E0 / e_mass
    gamma2 = gamma * gamma
    beta = sqrt(1.0 - 1.0/gamma2)
    U0 = Cgamma / 2.0 / pi * E0**4 * rp.i2
    Jx = 1.0 - rp.i4/rp.i2
    Jz = 1.0
    Je = 2.0 + rp.i4/rp.i2
    damping_partition_numbers = numpy.array([Jx, Jz, Je])
    ct = 2.0 * E0 / U0 * revolution_period
    rp.E0 = E0
    rp.U0 = U0
    emitx = Cq * gamma2 * rp.i5 / Jx / rp.i2
    rp.emittances = numpy.array([emitx, nan, nan])
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
    ringtunes, _ = numpy.modf(ring.periodicity * tunes)
    rp.fulltunes = ring.periodicity * twiss[-1].mu / 2.0 / pi
    rp.tunes = numpy.concatenate((ringtunes, (nus,)))
    rp.alphac = alphac
    rp.etac = etac
    return rp


# noinspection PyPep8Naming
def envelope_parameters(ring, params=None):
    """Compute ring parameters from ohmi_envelope

    INPUT
        ring            Lattice object.

    KEYWORD
        params=None     RingParam object to be updated.

    OUTPUT
        params          RingParam object. The computed attributes are,

            tunes6          (3,) fractional (H, V, Long.) tunes (6D motion)
            emittances      (3,) Mode emittances
            J               (3,) Damping partition numbers
            Tau             (3,) Damping times [s]
            sigma_e         Energy spread
            sigma_l         Bunch length [m]
            voltage         Total accelerating voltage [V]
            phi_s           Synchrotron phase [rad]
            f_s             Synchrotron frequency [Hz]
    """
    rp = RingParameters() if params is None else params
    emit0, beamdata, emit = ring.ohmi_envelope()
    voltage = ring.voltage
    rp.E0 = ring.energy
    rp.U0 = ring.energy_loss
    revolution_period = ring.circumference / c
    rp.Tau = revolution_period / beamdata.damping_rates / ring.periodicity
    alpha = 1.0 / rp.Tau
    rp.J = 4.0 * alpha / numpy.sum(alpha)
    rp.tunes6, _ = numpy.modf(ring.periodicity * beamdata.tunes)
    rp.phi_s = pi - asin(rp.U0 / voltage)
    rp.voltage = voltage
    rp.f_s = rp.tunes6[2] / revolution_period
    rp.emittances = beamdata.mode_emittances
    rp.sigma_e = sqrt(emit0.r66[4, 4])
    rp.sigma_l = sqrt(emit0.r66[5, 5])
    return rp


Lattice.radiation_parameters = radiation_parameters
Lattice.envelope_parameters = envelope_parameters
