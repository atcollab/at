from math import pi, sqrt, asin
import numpy
from numpy import nan
from scipy.constants import physical_constants as cst

_e_radius = cst['classical electron radius'][0]
_e_mass = cst['electron mass energy equivalent in MeV'][0]
_c = cst['speed of light in vacuum'][0]
_hbar = cst['Planck constant over 2 pi times c in MeV fm'][0]

Cgamma = 4.0E9 * pi * _e_radius / 3.0 / pow(_e_mass, 3)     # m / GeV^3
Cq = 55 / 32 / sqrt(3) * _hbar / _e_mass * 1.0e-15          # m


_fields = {
    'tunes':        '              Frac. tunes: {0}',
    'fulltunes':    '                    Tunes: {0}',
    'chromaticities': '           Chromaticities: {0}',
    'E0':           '                   Energy: {0:e} eV',
    'U0':           '       Energy loss / turn: {0:e} eV',
    'i1':           ' Radiation integrals - I1: {0} m',
    'i2':           '                       I2: {0} m^-1',
    'i3':           '                       I3: {0} m^-2',
    'i4':           '                       I4: {0} m^-1',
    'i5':           '                       I5: {0} m^-1',
    'emittances':   '          Mode emittances: {0}',
    'J':            'Damping partition numbers: {0}',
    'Tau':          '            Damping times: {0} s',
    'sigma_e':      '            Energy spread: {0:g}',
    'sigma_l':      '             Bunch length: {0:g} m',
    'voltage':      '         Cavities voltage: {0} V',
    'phi_s':        '        Synchrotron phase: {0:g} rd',
    'f_s':          '    Synchrotron frequency; {0:g} Hz'
}


class RingParameters(object):
    """Class for pretty printing the ring properties"""
    def __str__(self):
        vrs = vars(self).copy()
        vals = [(_fields[k], vrs.pop(k, None)) for k in _fields]
        # Predefined attributes
        lines = [f.format(v) for f, v in vals if v is not None]
        # Other attributes
        lines += ['{0:>25}: {1}'.format(k, getattr(self, k)) for k in vrs]
        return '\n'.join(lines)


# noinspection PyPep8Naming
def radiation_parameters(ring, dp=0.0, params=None):
    """Compute ring parameters from the radiation integrals. N

    INPUT
        ring            Lattice object.
        dp=0.0          momentum deviation

    KEYWORD
        params=None     RingParam object to be updated

    OUTPUT
        params          RingParam object. The computed attributes are,

            tunes           (3,) fractional (H, V, Long.) tunes
            fulltunes       (3,) full tunes
            chromaticities  (2,) H, V Chromaticities
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
    circumference = ring.get_s_pos(len(ring))[0] * ring.periodicity
    revolution_period = circumference / _c
    voltage = ring.voltage
    E0 = ring.energy
    gamma = 1.0e-6 * E0 / _e_mass
    U0 = 1.0e-27 * Cgamma / 2.0 / pi * E0**4 * rp.i2
    Jx = 1.0 - rp.i4/rp.i2
    Jz = 1.0
    Je = 2.0 + rp.i4/rp.i2
    damping_partition_numbers = numpy.array([Jx, Jz, Je])
    ct = 2.0 * E0 / U0 * circumference / _c
    rp.E0 = E0
    rp.U0 = U0
    emitx = Cq * gamma * gamma * rp.i5 / Jx / rp.i2
    rp.emittances = numpy.array([emitx, nan, nan])
    alphac = rp.i1 / circumference
    rp.phi_s = pi - asin(U0 / voltage)
    nus = sqrt(alphac * ring.harmonic_number *
               sqrt(voltage*voltage - U0*U0) / 2.0 / pi / E0)
    rp.voltage = voltage
    rp.f_s = nus / revolution_period
    rp.Tau = ct / damping_partition_numbers
    rp.J = damping_partition_numbers
    rp.sigma_e = sqrt(Cq * gamma * gamma * rp.i3 / Je / rp.i2)
    rp.sigma_l = alphac * circumference / 2.0 / pi / nus * rp.sigma_e
    ringtunes, _ = numpy.modf(ring.periodicity * tunes)
    rp.fulltunes = ring.periodicity * twiss[-1].mu / 2.0 / pi
    rp.tunes = numpy.concatenate((ringtunes, (nus,)))
    return rp


# noinspection PyPep8Naming
def envelope_parameters(ring, params=None):
    """Compute ring parameters from ohmi_envelope

    INPUT
        ring            Lattice object.

    KEYWORD
        params=None     RingParam object to be updated

    OUTPUT
        params          RingParam object. The computed attributes are,

            tunes           (3,) fractional (H, V, Long.) tunes
            E0              Energy [eV]
            U0              Energy loss / turn [eV]
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
    circumference = ring.get_s_pos(len(ring))[0] * ring.periodicity
    voltage = ring.voltage
    rp.E0 = ring.energy
    rp.U0 = ring.energy_loss
    revolution_period = circumference / _c
    rp.Tau = revolution_period / beamdata.damping_rates / ring.periodicity
    alpha = 1.0 / rp.Tau
    rp.J = 4.0 * alpha / numpy.sum(alpha)
    rp.tunes, _ = numpy.modf(ring.periodicity * beamdata.tunes)
    rp.phi_s = pi - asin(rp.U0 / voltage)
    rp.voltage = voltage
    rp.f_s = rp.tunes[2] / revolution_period
    rp.emittances = beamdata.mode_emittances
    rp.sigma_e = sqrt(emit0.r66[4, 4])
    rp.sigma_l = sqrt(emit0.r66[5, 5])
    return rp
