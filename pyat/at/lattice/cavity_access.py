from enum import Enum
import numpy
from .elements import RFCavity
from .utils import AtError, checktype, make_copy, get_cells
from .lattice_object import Lattice

__all__ = ['get_rf_frequency', 'get_rf_voltage', 'set_rf_voltage',
           'get_rf_timelag', 'set_rf_timelag', 'set_cavity', 'Frf']


class Frf(Enum):
    """Enum class for frequency setting"""
    NOMINAL = 'nominal'


def _select_cav(ring, cavpts):
    """Select the cavities"""
    if cavpts is None:
        try:
            cavpts = ring.cavpts
        except AttributeError:
            cavpts = get_cells(ring, checktype(RFCavity))
    return cavpts


def _singlev(values, attr):
    """Return the single value"""
    vals = numpy.unique(values)
    if len(vals) > 1:
        raise AtError('{0} not equal for all cavities'.format(attr))
    return vals[0]


# noinspection PyUnusedLocal
def _sumv(values, attr):
    """Return the sum of values"""
    return numpy.sum(values)


def _fundmask(ring, cavpts):
    freqs = numpy.array([cav.Frequency for cav in ring.select(cavpts)])
    if len(freqs) < 1:
        raise AtError('No cavity found in the lattice')
    mask = (freqs == min(freqs))
    return None if numpy.all(mask) else mask


def _get_cavity(ring, attr, fsingle, cavpts=None, array=False):
    cavpts = _select_cav(ring, cavpts)
    vcell = numpy.array([getattr(cav, attr) for cav in ring.select(cavpts)])
    if array:
        return vcell
    else:
        fundmask = _fundmask(ring, cavpts)
        if fundmask is not None:
            vcell = vcell[fundmask]
        return fsingle(vcell, attr)


def get_rf_frequency(ring, **kwargs):
    """Return the RF frequency
    KEYWORDS
        cavpts=None   Cavity location.
                      If None, look for ring.cavpts, or otherwise take all
                      cavities.
        array=False   If False, return the frequency of the selected cavities
                      with the lowest frequency.
                      If True, return the frequency of all selected cavities
    """
    return _get_cavity(ring, 'Frequency', _singlev, **kwargs)


def get_rf_voltage(ring, **kwargs):
    """Return the total RF voltage (full ring)
    KEYWORDS
        cavpts=None   Cavity location.
                      If None, look for ring.cavpts, or otherwise take all
                      cavities.
        array=False   If False, return the sum of the voltage of the selected
                      cavities with the lowest frequency.
                      If True, return the voltage of all the selected cavities.
    """
    vcell = _get_cavity(ring, 'Voltage', _sumv, **kwargs)
    return ring.periodicity * vcell


def get_rf_timelag(ring, **kwargs):
    """Return the RF time lag
    KEYWORDS
        cavpts=None   Cavity location.
                      If None, look for ring.cavpts, or otherwise take all
                      cavities.
        array=False   If False, return the time lag of the cavities with the
                      lowest frequency.
                      If True, return the time lag of all the selected cavities.
    """
    return _get_cavity(ring, 'TimeLag', _singlev, **kwargs)


def set_rf_voltage(ring, voltage, **kwargs):
    """Set the RF voltage for the full ring

    PARAMETERS
        ring            lattice description
        voltage         RF voltage [V]

    KEYWORDS
        cavpts=None     If None, look for ring.cavpts, or otherwise take all
                        cavities.
        array=False     If False, the voltages of all cavities are scaled to
                        reach the specified value on the selected cavities with
                        the lowest frequency.
                        If True, directly apply voltage to the selected
                        cavities. The value must be broadcastable to the number
                        of cavities.
        copy=False      If True, returns a shallow copy of ring with new
                        cavity elements. Otherwise, modify ring in-place
    """
    return set_cavity(ring, Voltage=voltage, **kwargs)


def set_rf_timelag(ring, timelag, **kwargs):
    """Set the RF time lag

    PARAMETERS
        ring            lattice description
        timelag         RF time shift (-ct) [m]

    KEYWORDS
        cavpts=None     If None, look for ring.cavpts, or otherwise take all
                        cavities.
        array=False     If False, timelag is applied to the selected cavities
                        with the lowest frequency. The timelag of all the
                        other selected cavities is shifted by the same amount.
                        If True, directly apply timelag to the selected
                        cavities. The value must be broadcastable to the number
                        of cavities.
        copy=False      If True, returns a shallow copy of ring with new
                        cavity elements. Otherwise, modify ring in-place
    """
    return set_cavity(ring, TimeLag=timelag, **kwargs)


# noinspection PyPep8Naming
def set_cavity(ring, Voltage=None, Frequency=None, TimeLag=None, cavpts=None,
               copy=False, array=False):
    """
    Set the parameters of the RF cavities

    PARAMETERS
        ring                lattice description

    KEYWORDS
        Frequency=None  RF frequency [Hz]
        Voltage=None    RF voltage [V]
        TimeLag=None    RF time shift [-ct]
        cavpts=None     Cavity location. If None, look for ring.cavpts, or
                        otherwise take all cavities
        array=False     If False, the value is applied as described for
                        set_rf_voltage, set_rf_timelag and set_rf_frequency
                        If True, directly apply the value to the selected
                        cavities. The value must be broadcastable to the number
                        of cavities.
        copy=False      If True, returns a shallow copy of ring with new
                        cavity elements. Otherwise, modify ring in-place
    """
    # noinspection PyShadowingNames
    def getv(ring, attr, refpts):
        values = numpy.array([getattr(c, attr) for c in ring.select(refpts)])
        valfund = values[fundmask]
        return values, valfund

    cavpts = _select_cav(ring, cavpts)
    fundmask = None if array else _fundmask(ring, cavpts)

    modif = {}
    if Frequency is not None:
        if Frequency is Frf.NOMINAL:
            Frequency = ring.revolution_frequency * ring.harmonic_number
        if fundmask is not None:
            vall, vmain = getv(ring, 'Frequency', cavpts)
            vmain = vmain[0]
            if vmain == 0.0:
                vall[fundmask] = Frequency
                Frequency = vall
            else:
                Frequency *= vall/vmain
        modif['Frequency'] = Frequency

    if TimeLag is not None:
        if fundmask is not None:
            vall, vmain = getv(ring, 'TimeLag', cavpts)
            TimeLag += vall-_singlev(vmain, 'TimeLag')
        modif['TimeLag'] = TimeLag
        
    if Voltage is not None:
        vcell = Voltage / ring.periodicity
        if fundmask is not None:
            vall, vmain = getv(ring, 'Voltage', cavpts)
            vmain = numpy.sum(vmain)
            if vmain == 0.0:
                vall[fundmask] = vcell/numpy.count_nonzero(fundmask)
                vcell = vall
            else:
                vcell *= vall/vmain
        modif['Voltage'] = vcell

    # noinspection PyShadowingNames
    @make_copy(copy)
    def apply(ring, cavpts, modif):
        ncavs = ring.refcount(cavpts)
        for attr in modif.keys():
            try:
                values = numpy.broadcast_to(modif[attr], (ncavs,))
            except ValueError:
                raise AtError('set_cavity args should be either scalar or '
                              'a ({0},) vector'.format(ncavs))
            for val, cavity in zip(values, ring.select(cavpts)):
                cavity.update({attr: val})

    return apply(ring, cavpts, modif)


Lattice.get_rf_voltage = get_rf_voltage
Lattice.get_rf_frequency = get_rf_frequency
Lattice.get_rf_timelag = get_rf_timelag
Lattice.set_rf_voltage = set_rf_voltage
Lattice.set_rf_timelag = set_rf_timelag
Lattice.set_cavity = set_cavity
Lattice.rf_voltage = property(get_rf_voltage, set_rf_voltage,
                              doc="RF voltage of the full ring [V]")
Lattice.rf_timelag = property(get_rf_timelag, set_rf_timelag,
                              doc="Time lag of the fundamental mode [m]")
