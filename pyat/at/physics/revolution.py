from ..lattice import Lattice, get_rf_frequency, check_6d, get_s_pos
from ..lattice import DConstant
from ..constants import clight
from ..tracking import lattice_pass
from .orbit import find_orbit4
import numpy
import functools
from typing import Optional

__all__ = ['frequency_control', 'get_mcf', 'get_slip_factor',
           'get_revolution_frequency', 'set_rf_frequency']


def frequency_control(func):
    r"""Function to be used as decorator for ``func(ring, *args, **kwargs)``

    If ``ring.is_6d`` is :py:obj:`True` **and** *dp*, *dct* or *df*
    is specified in *kwargs*, make a copy of *ring* with a modified
    RF frequency, remove *dp*, *dct* or *df* from *kwargs* and call
    *func* with the modified *ring*.

    If ``ring.is_6d`` is :py:obj:`False` **or** no *dp*, *dct* or
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
                ring = set_rf_frequency(ring, **momargs, copy=True)
        return func(ring, *args, **kwargs)
    return wrapper


@check_6d(False)
def get_mcf(ring: Lattice, dp: Optional[float] = 0.0,
            keep_lattice: bool = False, **kwargs) -> float:
    r"""Compute the momentum compaction factor :math:`\alpha`

    Parameters:
        ring:           Lattice description (``ring.is_6d`` must be
          :py:obj:`False`)
        dp:             Momentum deviation. Defaults to :py:obj:`None`
        keep_lattice:   Assume no lattice change since the previous tracking.
          Default: :py:obj:`False`

    Keyword Args:
        DPStep (Optional[float]):       Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`

    Returns:
        mcf (float):    Momentum compaction factor :math:`\alpha`
    """
    dp_step = kwargs.pop('DPStep', DConstant.DPStep)
    fp_a, _ = find_orbit4(ring, dp=dp - 0.5*dp_step, keep_lattice=keep_lattice)
    fp_b, _ = find_orbit4(ring, dp=dp + 0.5*dp_step, keep_lattice=True)
    fp = numpy.stack((fp_a, fp_b),
                     axis=0).T  # generate a Fortran contiguous array
    b = numpy.squeeze(lattice_pass(ring, fp, keep_lattice=True), axis=(2, 3))
    ring_length = get_s_pos(ring, len(ring))
    alphac = (b[5, 1] - b[5, 0]) / dp_step / ring_length[0]
    return alphac


def get_slip_factor(ring: Lattice, **kwargs) -> float:
    r"""Compute the slip factor :math:`\eta`

    Parameters:
        ring:           Lattice description (``ring.is_6d`` must be
          :py:obj:`False`)

    Keyword Args:
        dp (float):     Momentum deviation
        DPStep (float): Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`

    Returns:
        eta (float):    Slip factor :math:`\eta`
    """
    gamma = ring.gamma
    etac = (1.0/gamma/gamma - get_mcf(ring, **kwargs))
    return etac


def get_revolution_frequency(ring: Lattice,
                             dp: float = None,
                             dct: float = None,
                             df: float = None) -> float:
    """Compute the revolution frequency of the full ring [Hz]

    Parameters:
        ring:       Lattice description
        dp:         Momentum deviation. Defaults to :py:obj:`None`
        dct:        Path lengthening. Defaults to :py:obj:`None`
        df:         Deviation of RF frequency. Defaults to :py:obj:`None`

    Returns:
        frev:       Revolution frequency [Hz]
    """
    lcell = ring.cell_length
    cell_frev = ring.beta * clight / lcell
    if dct is not None:
        cell_frev *= lcell / (lcell + dct)
    elif dp is not None:
        # Find the path lengthening for dp
        rnorad = ring.disable_6d(copy=True) if ring.is_6d else ring
        orbit = lattice_pass(rnorad, rnorad.find_orbit4(dp=dp)[0])
        dct = numpy.squeeze(orbit)[5]
        cell_frev *= lcell / (lcell + dct)
    elif df is not None:
        cell_frev += df / ring.cell_harmnumber
    return cell_frev / ring.periodicity


def set_rf_frequency(ring: Lattice, frequency: float = None,
                     dp: float = None, dct: float = None, df: float = None,
                     **kwargs):
    """Set the RF frequency

    Parameters:
        ring:           Lattice description
        frequency:      RF frequency [Hz]. Default: nominal frequency.
        dp:             Momentum deviation. Defaults to :py:obj:`None`
        dct:            Path lengthening. Defaults to :py:obj:`None`
        df:             Deviation of RF frequency. Defaults to :py:obj:`None`

    Keyword Args:
        cavpts (Optional[Refpts]):  If :py:obj:`None`, look for ring.cavpts, or
          otherwise take all cavities.
        array (Optional[bool]):     If :py:obj:`False` (default), *frequency*
          is applied to the selected cavities with the lowest frequency. The
          frequency of all the other selected cavities is scaled by the same
          ratio.

          If :py:obj:`True`, directly apply *frequency* to the selected
          cavities. The value must be broadcastable to the number of cavities.
        copy (Optional[bool]):     If :py:obj:`True`, returns a shallow copy of
          *ring* with new cavity elements. Otherwise (default), modify
          *ring* in-place
    """
    if frequency is None:
        frequency = ring.get_revolution_frequency(dp=dp, dct=dct, df=df) \
                    * ring.harmonic_number
    return ring.set_cavity(Frequency=frequency, **kwargs)


Lattice.mcf = property(get_mcf, doc="Momentum compaction factor")
Lattice.slip_factor = property(get_slip_factor, doc="Slip factor")
Lattice.get_revolution_frequency = get_revolution_frequency
Lattice.get_mcf = get_mcf
Lattice.get_slip_factor = get_slip_factor
Lattice.set_rf_frequency = set_rf_frequency
Lattice.rf_frequency = property(get_rf_frequency, set_rf_frequency,
    doc="Fundamental RF frequency [Hz]. The special value "
        ":py:class:`Frf.NOMINAL <.Frf>` means nominal frequency.")
