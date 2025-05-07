"""Access to properties of RF cavities

See :doc:`/p/howto/CavityControl` for a guide on using these functions.
"""

from __future__ import annotations

import functools
from enum import Enum

import numpy as np

from .elements import RFCavity
from .lattice_object import Lattice
from .utils import AtError, make_copy, Refpts

__all__ = [
    "frequency_control",
    "get_rf_frequency",
    "set_rf_frequency",
    "get_rf_voltage",
    "set_rf_voltage",
    "get_rf_timelag",
    "set_rf_timelag",
    "set_cavity",
    "Frf",
]


def frequency_control(func):
    r"""Function to be used as decorator for
    :pycode:`func(ring, *args, **kwargs)`

    If :pycode:`ring.is_6d` is :py:obj:`True` **and** *dp*, *dct* or *df*
    is specified in *kwargs*, make a copy of *ring* with a modified
    RF frequency, remove *dp*, *dct* or *df* from *kwargs* and call
    *func* with the modified *ring*.

    If :pycode:`ring.is_6d` is :py:obj:`False` **or** no *dp*, *dct* or
    *df* is specified in *kwargs*, *func* is called with *ring* unchanged.

    Example:

        >>> @frequency_control
        ... def func(ring, *args, dp=None, dct=None, df=None, **kwargs):
        ...     pass
    """

    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs):
        if ring.is_6d:
            momargs = {}
            for key in ["dp", "dct", "df"]:
                v = kwargs.pop(key, None)
                if v is not None:
                    momargs[key] = v
            if len(momargs) > 0:
                ring = set_rf_frequency(ring, copy=True, **momargs)
        return func(ring, *args, **kwargs)

    return wrapper


class Frf(Enum):
    """Enum class for frequency setting"""

    #:  Constant used as a frequency value, standing for the nominal frequency
    NOMINAL = "nominal"


def _selected_cavities(ring: Lattice, cavpts: Refpts) -> list[RFCavity]:
    """Select the cavities"""
    if cavpts is None:
        try:
            cavpts = ring.cavpts
        except AttributeError:
            cavpts = RFCavity
    return cavpts


def _singlev(values, attr):
    """Return the single value"""
    vals = np.unique(values)
    if len(vals) > 1:
        raise AtError(f"{attr} not equal for all cavities")
    return vals[0]


# noinspection PyUnusedLocal
def _sumv(values, attr):
    """Return the sum of values"""
    return np.sum(values)


def _fundmask(cavities: list[RFCavity]):
    freqs = np.array([cav.Frequency for cav in cavities])
    if len(freqs) < 1:
        raise AtError("No cavity found in the lattice")
    mask = freqs == np.min(freqs)
    return mask


def _get_cavity(ring: Lattice, attr, fsingle, cavpts=None, array=False):
    cavpts = _selected_cavities(ring, cavpts)
    cavities = list(ring.select(cavpts))
    vcell = np.array([getattr(cav, attr) for cav in cavities])
    if array:
        return vcell
    else:
        fundmask = _fundmask(cavities)
        if fundmask is not None:
            vcell = vcell[fundmask]
        return fsingle(vcell, attr)


def get_rf_frequency(ring: Lattice, **kwargs) -> float:
    """Return the RF frequency

    Parameters:
        ring:               Lattice description

    Keyword Arguments:
        cavpts (Refpts | None):  If :py:obj:`None`, look for :pycode:`ring.cavpts`, or
          otherwise take all cavities. The main cavities are those with the lowest
          frequency among the selected ones.
        array (bool):       If :py:obj:`False`, return the frequency of the main
          cavities (selected cavities with the lowest frequency).

          If :py:obj:`True`, return the frequency of each selected cavity.

    Returns:
        rf_freq:    RF frequency [Hz]
    """
    return _get_cavity(ring, "Frequency", _singlev, **kwargs)


def get_rf_voltage(ring: Lattice, **kwargs) -> float:
    """Return the total RF voltage (full ring including periodicity)

    Parameters:
        ring:               Lattice description

    Keyword Arguments:
        cavpts (Refpts | None):  If :py:obj:`None`, look for :pycode:`ring.cavpts`, or
          otherwise take all cavities. The main cavities are those with the lowest
          frequency among the selected ones.
        array (bool):       If :py:obj:`False`, return the sum of the voltage of the
          main cavities (selected cavities with the lowest frequency).

          If :py:obj:`True`, return the voltage of each selected cavity.

    Returns:
        rf_v:       Total voltage (full ring including periodicity) [V]
    """
    vcell = _get_cavity(ring, "Voltage", _sumv, **kwargs)
    return ring.periodicity * vcell


def get_rf_timelag(ring: Lattice, **kwargs) -> float:
    """Return the RF time lag

    Parameters:
        ring:           Lattice description

    Keyword Arguments:
        cavpts (Refpts | None):  If :py:obj:`None`, look for :pycode:`ring.cavpts`, or
          otherwise take all cavities. The main cavities are those with the lowest
          frequency among the selected ones.
        array (bool):       If :py:obj:`False`, return the timelag of the main cavities
          (selected cavities with the lowest frequency).

          If :py:obj:`True`, return the timelag of each selected cavity.

    Returns:
        rf_timelag: RF time lag as path lengthening [m]
    """
    return _get_cavity(ring, "TimeLag", _singlev, **kwargs)


def set_rf_frequency(
    ring: Lattice, frequency: float | Frf = Frf.NOMINAL, **kwargs
) -> Lattice | None:
    """Set the RF frequency

    Parameters:
        ring:           Lattice description
        frequency:      RF frequency [Hz]. Use :py:obj:`.Frf.NOMINAL` to specify the
          nominal frequency corresponding to the given off-momentum.

    Keyword Args:
        dp (float | None):  Momentum deviation. Defaults to :py:obj:`None`
        dct (float | None): Path lengthening. Defaults to :py:obj:`None`
        df (float | None):  Deviation of RF frequency. Defaults to :py:obj:`None`
        cavpts (Refpts | None):  If :py:obj:`None`, look for :pycode:`ring.cavpts`, or
          otherwise take all cavities.
        array (bool):       If :py:obj:`False` (default), *frequency* is applied to
          the selected cavities with the lowest frequency. The frequency of all the
          other selected cavities is scaled by the same ratio.

          If :py:obj:`True`, directly apply *frequency* to the selected cavities. The
          value must be broadcastable to the number of cavities.
        copy (bool):        If :py:obj:`True`, returns a shallow copy of *ring* with
          new cavity elements. Otherwise (default), modify *ring* in-place.
    """
    return set_cavity(ring, Frequency=frequency, **kwargs)


def set_rf_voltage(ring: Lattice, voltage: float, **kwargs) -> Lattice | None:
    """Set the RF voltage for the full ring

    Parameters:
        ring:           Lattice description
        voltage:        Total RF voltage (full ring taking periodicity into account) [V]

    Keyword Arguments:
        cavpts (Refpts | None):  If :py:obj:`None`, look for :pycode:`ring.cavpts`, or
          otherwise take all cavities.
        array (bool):       If :py:obj:`False` (default), the voltage of the main
          cavities is scaled to achieve the desired *voltage*. The voltage of the other
          cavities is scaled by the same ratio.

          If :py:obj:`True`, directly apply *voltage* to the selected cavities. The
          value must be broadcastable to the number of cavities.
        copy (bool):        If :py:obj:`True`, returns a shallow copy of *ring* with
          new cavity elements. Otherwise (default), modify *ring* in-place.
    """
    return set_cavity(ring, Voltage=voltage, **kwargs)


def set_rf_timelag(ring: Lattice, timelag: float, **kwargs) -> Lattice | None:
    """Set the RF time lag

    Parameters:
        ring:           Lattice description
        timelag:        RF time lag, expressed as path lengthening [m]

    Keyword Arguments:
        cavpts (Refpts | None):  If :py:obj:`None`, look for :pycode:`ring.cavpts`, or
          otherwise take all cavities.
        array (bool):       If :py:obj:`False` (default), *timelag* is applied to the
          main cavities. The timelag of the other cavities is shifted by the same
          amount.

          If :py:obj:`True`, directly apply *timelag* to the selected cavities. The
          value must be broadcastable to the number of cavities.
        copy (bool):    If :py:obj:`True`, returns a shallow copy of *ring* with new
          cavity elements. Otherwise (default), modify *ring* in-place.
    """
    return set_cavity(ring, TimeLag=timelag, **kwargs)


# noinspection PyPep8Naming
def set_cavity(
    ring: Lattice,
    Voltage: float | None = None,
    Frequency: float | Frf | None = None,
    TimeLag: float | None = None,
    cavpts: Refpts | None = None,
    copy: bool = False,
    array: bool = False,
    **kwargs,
) -> Lattice | None:
    """
    Set the parameters of the RF cavities

    Parameters:
        ring:           lattice description
        Frequency:      RF frequency [Hz]. Use :py:obj:`.Frf.NOMINAL` to specify the
          nominal frequency corresponding to the given off-momentum.
        Voltage:        Total RF voltage [V]
        TimeLag:        RF time shift expressed as path lengthening [m]
        cavpts:         Cavity location. If None, look for :pycode:`ring.cavpts`, or
          otherwise take all cavities
        array:          If :py:obj:`False`, the value is applied as described for
          :py:func:`set_rf_voltage`, :py:func:`set_rf_timelag` and
          :py:func:`set_rf_frequency`.

          If :py:obj:`True`, directly apply the values to the selected cavities. The
          values must be broadcastable to the number of cavities.
        copy (bool):    If :py:obj:`True`, returns a shallow copy of *ring* with new
          cavity elements. Otherwise (default), modify *ring* in-place.
    """

    # noinspection PyShadowingNames
    def getv(cavities: list[RFCavity], attr: str):
        values = np.array([getattr(el, attr) for el in cavities])
        valfund = values[fundmask]
        return values, valfund

    def nominal(elem: RFCavity, cell_frev: float) -> float:
        cell_h = getattr(elem, "HarmNumber", 0)
        if cell_h == 0:
            cell_h = round(elem.Frequency / cell_frev)
        return cell_h * cell_frev

    # noinspection PyShadowingNames
    @make_copy(copy)
    def apply(rg: Lattice, cavpts, modif):
        cavities = list(rg.select(cavpts))
        ncavs = len(cavities)
        for attr, values in modif.items():
            try:
                values = np.broadcast_to(values, (ncavs,))
            except ValueError:
                raise ValueError(
                    f"{attr} should be either scalar or a ({ncavs},) array"
                ) from None
            for val, cavity in zip(values, cavities):
                setattr(cavity, attr, val)

    cavpts = _selected_cavities(ring, cavpts)
    cavities = list(ring.select(cavpts))
    fundmask = _fundmask(cavities)

    modif = {}
    if Frequency is not None:
        if Frequency is Frf.NOMINAL:
            cell_f = ring.get_revolution_frequency(**kwargs) * ring.periodicity
            Frequency = np.array([nominal(el, cell_f) for el in cavities])
        elif not array:
            vall, vmain = getv(cavities, "Frequency")
            vmain = vmain[0]
            if vmain == 0.0:
                vall[fundmask] = Frequency
                Frequency = vall
            else:
                Frequency *= vall / vmain
        modif["Frequency"] = Frequency

    if TimeLag is not None:
        if not array:
            vall, vmain = getv(cavities, "TimeLag")
            TimeLag += vall - _singlev(vmain, "TimeLag")
        modif["TimeLag"] = TimeLag

    if Voltage is not None:
        vcell = Voltage / ring.periodicity
        if not array:
            vall, vmain = getv(cavities, "Voltage")
            vmain = np.sum(vmain)
            if vmain == 0.0:
                vall[fundmask] = vcell / np.count_nonzero(fundmask)
                vcell = vall
            else:
                vcell *= vall / vmain
        modif["Voltage"] = vcell

    return apply(ring, cavpts, modif)


Lattice.get_rf_voltage = get_rf_voltage
Lattice.get_rf_frequency = get_rf_frequency
Lattice.get_rf_timelag = get_rf_timelag
Lattice.set_rf_voltage = set_rf_voltage
Lattice.set_rf_frequency = set_rf_frequency
Lattice.set_rf_timelag = set_rf_timelag
Lattice.set_cavity = set_cavity
Lattice.rf_voltage = property(
    get_rf_voltage, set_rf_voltage, doc="RF voltage of the full ring [V]"
)
Lattice.rf_timelag = property(
    get_rf_timelag, set_rf_timelag, doc="Time lag of the fundamental mode [m]"
)
Lattice.rf_frequency = property(
    get_rf_frequency,
    set_rf_frequency,
    doc="Fundamental RF frequency [Hz]. The special value "
    ":py:class:`Frf.NOMINAL <.Frf>` means nominal frequency.",
)
