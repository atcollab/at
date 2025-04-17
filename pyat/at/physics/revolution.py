"""Revolution frequency, momentum compaction factor, slip factor"""
from __future__ import annotations

__all__ = ["get_mcf", "get_slip_factor", "get_revolution_frequency"]

import numpy as np

from .orbit import find_orbit4
from ..constants import clight
from ..lattice import DConstant, Lattice, check_6d, get_s_pos
from ..tracking import internal_lpass


@check_6d(False)
def get_mcf(
    ring: Lattice,
    dp: float = 0.0,
    keep_lattice: bool = False,
    fit_order: int = 1,
    n_step: int = 2,
    **kwargs,
) -> float:
    r"""Compute the momentum compaction factor :math:`\alpha`

    Parameters:
        ring:           Lattice description (:pycode:`ring.is_6d` must be
          :py:obj:`False`)
        dp:             Momentum deviation
        keep_lattice:   Assume no lattice change since the previous tracking.
        fit_order:      Maximum momentum compaction factor order to be fitted.
          Default to 1, corresponding to the first-order momentum compaction factor.
        n_step:         Number of different calculated momentum deviations to be fitted
          with a polynomial.
          Default to 2.

    Keyword Args:
        DPStep (float): Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`

    Returns:
        mcf (float/array):    Momentum compaction factor :math:`\alpha`
          up to the order *fit_order*. Returns a float if *fit_order==1* otherwise
          returns an array.
    """
    if n_step < 2 * fit_order:
        raise ValueError(
            "Low number of steps, it is advised to have n_step >= 2*fit_order"
            + " for a better fit."
        )
    dp_step = kwargs.pop("DPStep", DConstant.DPStep)
    dp_samples = np.linspace(-dp_step / 2, dp_step / 2, n_step)
    fp_i = tuple(
        find_orbit4(ring, dp=dp + dp_i, keep_lattice=keep_lattice)[0]
        for dp_i in dp_samples
    )
    fp = np.stack(fp_i, axis=0).T  # generate a Fortran contiguous array
    b = np.squeeze(internal_lpass(ring, fp, keep_lattice=True), axis=(2, 3))
    ring_length = get_s_pos(ring, len(ring))
    p = np.polynomial.Polynomial.fit(b[4], b[5], deg=fit_order).convert().coef
    alphac = p[1:] / ring_length[0]
    return alphac[0] if len(alphac) < 2 else alphac


def get_slip_factor(ring: Lattice, **kwargs) -> float:
    r"""Compute the slip factor :math:`\eta=1/\gamma^2-\alpha`

    Parameters:
        ring:           Lattice description (:pycode:`ring.is_6d` must be
          :py:obj:`False`)

    Keyword Args:
        dp (float):     Momentum deviation
        DPStep (float): Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`

    Returns:
        eta (float):    Slip factor :math:`\eta`
    """
    gamma = ring.gamma
    etac = 1.0 / gamma / gamma - get_mcf(ring, **kwargs)
    return etac


def get_revolution_frequency(
    ring: Lattice, dp: float = None, dct: float = None, df: float = None
) -> float:
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
        orbit = internal_lpass(rnorad, rnorad.find_orbit4(dp=dp)[0])
        dct = np.squeeze(orbit)[5]
        cell_frev *= lcell / (lcell + dct)
    elif df is not None:
        cell_frev += df / ring.cell_harmnumber
    return cell_frev / ring.periodicity


Lattice.mcf = property(get_mcf, doc="Momentum compaction factor")
Lattice.slip_factor = property(get_slip_factor, doc="Slip factor")
Lattice.get_mcf = get_mcf
Lattice.get_slip_factor = get_slip_factor
Lattice.get_revolution_frequency = get_revolution_frequency
