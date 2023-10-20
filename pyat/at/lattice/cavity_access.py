"""Access to properties of RF cavities"""
from enum import Enum
import numpy
from typing import Optional
import functools
from .elements import RFCavity
from .utils import AtError, make_copy, Refpts, get_bool_index
from .lattice_object import Lattice

__all__ = ['frequency_control',
           'get_rf_frequency', 'get_rf_voltage', 'set_rf_voltage',
           'get_rf_timelag', 'set_rf_timelag', 'set_cavity', 'Frf']


def frequency_control(func):
    r""" Function to be used as decorator for
    :pycode:`func(ring, *args, **kwargs)`

    If :pycode:`ring.is_6d` is :py:obj:`True` **and** *dp*, *dct* or *df*
    is specified in *kwargs*, make a copy of *ring* with a modified
    RF frequency, remove *dp*, *dct* or *df* from *kwargs* and call
    *func* with the modified *ring*.

    If :pycode:`ring.is_6d` is :py:obj:`False` **or** no *dp*, *dct* or
    *df* is specified in *kwargs*, *func* is called unchanged.

    Examples:

        .. code-block:: python

            @frequency_control
            def func(ring, *args, dp=None, dct=None, **kwargs):
                pass
    """

    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs):
        if ring.is_6d:
            momargs = {}
            for key in ['dp', 'dct', 'df']:
                v = kwargs.pop(key, None)
                if v is not None:
                    momargs[key] = v
            if len(momargs) > 0:
                frequency = ring.get_revolution_frequency(**momargs) \
                            * ring.harmonic_number
                ring = set_cavity(ring, Frequency=frequency, copy=True)
        return func(ring, *args, **kwargs)

    return wrapper


class Frf(Enum):
    """Enum class for frequency setting"""
    NOMINAL = 'nominal'


def _select_cav(ring: Lattice, cavpts):
    """Select the cavities"""
    if cavpts is None:
        try:
            cavpts = ring.cavpts
        except AttributeError:
            cavpts = get_bool_index(ring, RFCavity)
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


def _fundmask(ring: Lattice, cavpts):
    freqs = numpy.array([cav.Frequency for cav in ring.select(cavpts)])
    if len(freqs) < 1:
        raise AtError('No cavity found in the lattice')
    mask = (freqs == min(freqs))
    return mask


def _get_cavity(ring: Lattice, attr, fsingle, cavpts=None, array=False):
    cavpts = _select_cav(ring, cavpts)
    vcell = numpy.array([getattr(cav, attr) for cav in ring.select(cavpts)])
    if array:
        return vcell
    else:
        fundmask = _fundmask(ring, cavpts)
        if fundmask is not None:
            vcell = vcell[fundmask]
        return fsingle(vcell, attr)


def get_rf_frequency(ring: Lattice, **kwargs) -> float:
    """Return the RF frequency

    Parameters:
        ring:           Lattice description

    Keyword Arguments:
        cavpts=None:    Cavity location. If None, look for ring.cavpts,
                          otherwise take all cavities.
        array=False:    If False, return the frequency of the selected cavities
                          with the lowest frequency.

                          If True, return the frequency of all selected cavities

    Returns:
        rf_freq:    RF frequency
    """
    return _get_cavity(ring, 'Frequency', _singlev, **kwargs)


def get_rf_voltage(ring: Lattice, **kwargs) -> float:
    """Return the total RF voltage (full ring)

    Parameters:
        ring:           Lattice description

    Keyword Arguments:
        cavpts=None:    Cavity location. If None, look for ring.cavpts,
                          otherwise take all cavities.
        array=False:    If False, return the frequency of the selected cavities
                          with the lowest frequency.

                          If True, return the frequency of all selected cavities

    Returns:
        rf_v:       Total RF voltage (full ring)
    """
    vcell = _get_cavity(ring, 'Voltage', _sumv, **kwargs)
    return ring.periodicity * vcell


def get_rf_timelag(ring: Lattice, **kwargs) -> float:
    """Return the RF time lag

    Parameters:
        ring:           Lattice description

    Keyword Arguments:
        cavpts=None:    Cavity location. If None, look for ring.cavpts,
                          otherwise take all cavities.
        array=False:    If False, return the frequency of the selected cavities
                          with the lowest frequency.

                          If True, return the frequency of all selected cavities

    Returns:
        rf_timelag: RF time lag
    """
    return _get_cavity(ring, 'TimeLag', _singlev, **kwargs)


def set_rf_voltage(ring: Lattice, voltage: float, **kwargs):
    """Set the RF voltage for the full ring

    Parameters:
        ring:           Lattice description
        voltage:        Total RF voltage (full ring)

    Keyword Arguments:
        cavpts=None:    Cavity location. If None, look for ring.cavpts,
                          otherwise take all cavities.
        array=False:    If False, return the frequency of the selected cavities
                          with the lowest frequency.

                          If True, return the frequency of all selected cavities
        copy=False:     If True, returns a shallow copy of ``ring`` with new
                          cavity elements. Otherwise, modify ``ring`` in-place.
    """
    return set_cavity(ring, Voltage=voltage, **kwargs)


def set_rf_timelag(ring: Lattice, timelag: float, **kwargs):
    """Set the RF time lag

    Parameters:
        ring:           Lattice description
        timelag:        RF time lag

    Keyword Arguments:
        cavpts=None:    Cavity location. If None, look for ring.cavpts,
                          otherwise take all cavities.
        array=False:    If False, return the frequency of the selected cavities
                          with the lowest frequency.

                          If True, return the frequency of all selected cavities
        copy=False:     If True, returns a shallow copy of ``ring`` with new
                          cavity elements. Otherwise, modify ``ring`` in-place.
    """
    return set_cavity(ring, TimeLag=timelag, **kwargs)


# noinspection PyPep8Naming
def set_cavity(ring: Lattice, Voltage: Optional[float] = None,
               Frequency: Optional[float] = None,
               TimeLag: Optional[float] = None,
               cavpts: Optional[Refpts] = None,
               copy: Optional[bool] = False,
               array: Optional[bool] = False):
    """
    Set the parameters of the RF cavities

    Parameters:
        ring:       lattice description
        Frequency:  RF frequency [Hz]
        Voltage:    RF voltage [V]
        TimeLag:    RF time shift [-ct]
        cavpts:     Cavity location. If None, look for ring.cavpts, or
                      otherwise take all cavities
        array:      If False, the value is applied as described for
                      set_rf_voltage, set_rf_timelag and set_rf_frequency
                      If True, directly apply the value to the selected
                      cavities. The value must be broadcastable to the number
                      of cavities.
        copy:       If True, returns a shallow copy of ring with new
                      cavity elements. Otherwise, modify ring in-place
    """
    # noinspection PyShadowingNames
    def getv(ring: Lattice, attr, refpts):
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
    def apply(ring: Lattice, cavpts, modif):
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
