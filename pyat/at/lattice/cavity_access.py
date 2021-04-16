from enum import Enum
import numpy
from scipy.constants import c as clight
from at.lattice import elements, AtError, checktype, make_copy, get_cells
from at.lattice import Lattice

__all__ = ['get_rf_frequency', 'get_rf_harmonic_number', 'get_rf_voltage',
           'set_rf_frequency', 'set_rf_harmonic_number', 'set_rf_voltage',
           'get_rf_timelag', 'set_rf_timelag', 'set_cavity', 'Frf']


class Frf(Enum):
    """Enum class for frequency setting"""
    NOMINAL = 'nominal'


def _get_rf_attr(ring, attr, cavpts=None):
    if cavpts is None:
        cavpts = checktype(elements.RFCavity)
    cavities = ring.select(cavpts)
    return numpy.array([getattr(cavity, attr) for cavity in cavities])


def _get_rf_unique_attr(ring, attr, cavpts=None):
    freq = numpy.unique(_get_rf_attr(ring, attr, cavpts=cavpts))
    if len(freq) == 0:
        raise AtError('No cavity found in the lattice')
    elif len(freq) > 1:
        raise AtError('{0} not equal for all cavities'.format(attr))
    else:
        return freq[0]


def get_rf_frequency(ring, cavpts=None):
    """Return the RF frequency
    KEYWORDS
        cavpts=None     location of RF cavities. If None, take all cavities.
                        This allows to ignore harmonic cavities
    """
    return _get_rf_unique_attr(ring, 'Frequency', cavpts=cavpts)


def get_rf_harmonic_number(ring, cavpts=None):
    """Return the RF harmonic number for the full ring
    KEYWORDS
        cavpts=None     location of RF cavities. If None, take all cavities.
                        This allows to ignore harmonic cavities
    """
    hcell = _get_rf_unique_attr(ring, 'HarmNumber', cavpts=cavpts)
    return ring.periodicity * hcell


def get_rf_voltage(ring, cavpts=None):
    """Return the total RF voltage for the full ring
    KEYWORDS
        cavpts=None     location of RF cavities. If None, take all cavities.
                        This allows to ignore harmonic cavities
    """
    _ = get_rf_frequency(ring, cavpts=cavpts)   # check single frequency
    vcell = sum(_get_rf_attr(ring, 'Voltage', cavpts=cavpts))
    return ring.periodicity * vcell


def get_rf_timelag(ring, cavpts=None):
    """Return the RF frequency
    KEYWORDS
        cavpts=None     location of RF cavities. If None, take all cavities.
                        This allows to ignore harmonic cavities
    """
    return _get_rf_unique_attr(ring, 'TimeLag', cavpts=cavpts)


def set_rf_frequency(ring, frequency, cavpts=None, copy=False):
    """Set the RF frequency

    PARAMETERS
        ring                lattice description
        frequency           RF frequency. The special enum value Frf.NOMINAL
                            sets the frequency to the nominal value, given
                            ring length and harmonic number.

    KEYWORDS
        cavpts=None         Cavity location. If None, use all cavities.
                            This allows to ignore harmonic cavities.
        copy=False          If True, returns a shallow copy of ring with new
                            cavity elements. Otherwise, modify ring in-place
    """
    return set_cavity(ring, Frequency=frequency, cavpts=cavpts, copy=copy)


def set_rf_harmonic_number(ring, harm_number, cavpts=None, copy=False):
    """Set the RF harmonic number

    PARAMETERS
        ring                lattice description
        harm_number         Harmonic number for the full ring

    KEYWORDS
        cavpts=None         Cavity location. If None, use all cavities.
                            This allows to ignore harmonic cavities.
        copy=False          If True, returns a shallow copy of ring with new
                            cavity elements. Otherwise, modify ring in-place
    """
    return set_cavity(ring, HarmNumber=harm_number, cavpts=cavpts, copy=copy)


def set_rf_voltage(ring, voltage, cavpts=None, copy=False):
    """Set the RF voltage

    PARAMETERS
        ring                lattice description
        voltage             Total RF voltage for the full ring

    KEYWORDS
        cavpts=None         Cavity location. If None, use all cavities.
                            This allows to ignore harmonic cavities.
        copy=False          If True, returns a shallow copy of ring with new
                            cavity elements. Otherwise, modify ring in-place
    """
    return set_cavity(ring, Voltage=voltage, cavpts=cavpts, copy=copy)


def set_rf_timelag(ring, timelag, cavpts=None, copy=False):
    """Set the RF time lag

    PARAMETERS
        ring                lattice description
        timelag

    KEYWORDS
        cavpts=None         Cavity location. If None, use all cavities.
                            This allows to ignore harmonic cavities.
        copy=False          If True, returns a shallow copy of ring with new
                            cavity elements. Otherwise, modify ring in-place
    """
    return set_cavity(ring, TimeLag=timelag, cavpts=cavpts, copy=copy)


# noinspection PyPep8Naming
def set_cavity(ring, Voltage=None, Frequency=None, HarmNumber=None,
               TimeLag=None, cavpts=None, copy=False):
    """
    PARAMETERS
        ring                lattice description

    KEYWORDS
        Frequency=None      RF frequency. The special enum value Frf.NOMINAL
                            sets the frequency to the nominal value, given
                            ring length and harmonic number.
        Voltage=None        RF voltage for the full ring.
        HarmNumber=None     Harmonic number for the full ring.
        TimeLag=None
        cavpts=None         Cavity location. If None, use all cavities.
                            This allows to ignore harmonic cavities.
        copy=False          If True, returns a shallow copy of ring with new
                            cavity elements. Otherwise, modify ring in-place
    """
    def get_harm(hcell):
        if hcell is not None:
            return hcell*n_periods
        else:
            return ring.get_rf_harmonic_number(cavpts=cavpts)

    if cavpts is None:
        cavpts = get_cells(ring, checktype(elements.RFCavity))
    n_cavities = ring.refcount(cavpts)
    if n_cavities < 1:
        raise AtError('No cavity found in the lattice')
    n_periods = ring.periodicity

    modif = {}
    if Frequency is not None:
        if Frequency is Frf.NOMINAL:
            Frequency = (clight / ring.circumference) * get_harm(HarmNumber)
        modif['Frequency'] = Frequency
    if HarmNumber is not None:
        modif['HarmNumber'] = HarmNumber / n_periods
    if TimeLag is not None:
        modif['TimeLag'] = TimeLag
    if Voltage is not None:
        modif['Voltage'] = Voltage / n_periods / n_cavities

    # noinspection PyShadowingNames
    @make_copy(copy)
    def apply(ring, cavpts, modif):
        for cavity in ring.select(cavpts):
            cavity.update(modif)

    return apply(ring, cavpts, modif)


Lattice.get_rf_voltage = get_rf_voltage
Lattice.get_rf_frequency = get_rf_frequency
Lattice.get_rf_harmonic_number = get_rf_harmonic_number
Lattice.get_rf_timelag = get_rf_timelag
Lattice.set_rf_voltage = set_rf_voltage
Lattice.set_rf_frequency = set_rf_frequency
Lattice.set_rf_harmonic_number = set_rf_harmonic_number
Lattice.set_rf_timelag = set_rf_timelag
Lattice.set_cavity = set_cavity
Lattice.voltage = property(get_rf_voltage, set_rf_voltage,
                           doc="Total RF voltage of the full ring [V]")
Lattice.harmonic_number = property(get_rf_harmonic_number,
                                   set_rf_harmonic_number,
                                   doc="Harmonic number of the full ring")
Lattice.frequency = property(get_rf_frequency, set_rf_frequency,
                             doc="RF frequency [Hz]")
