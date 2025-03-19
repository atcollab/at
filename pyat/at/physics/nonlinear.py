"""
Non-linear optics
"""

import numpy as np
from typing import Optional, Sequence
from scipy.special import factorial
from ..lattice import Element, Lattice
from ..tracking import internal_lpass
from .orbit import Orbit, find_orbit
from .linear import get_tune, get_chrom, linopt6
from .harmonic_analysis import get_tunes_harmonic

__all__ = ["detuning", "chromaticity", "gen_detuning_elem", "tunes_vs_amp"]


def tunes_vs_amp(
    ring: Lattice,
    amp: Optional[Sequence[float]] = None,
    dim: Optional[int] = 0,
    nturns: Optional[int] = 512,
    **kwargs,
) -> np.ndarray:
    r"""Generates a range of tunes for varying x, y, or z amplitudes

    Parameters:
        ring:               Lattice description
        amp:                Oscillation amplitude. If amp is None the
                            tunes are calculated from :py:func:`.linopt6`
                            else tracking is used
        dim:                Coordinate index of the oscillation
        nturns:             Number of turns
        method:             ``'laskar'`` or ``'fft'``. Default: ``'laskar'``
        num_harmonics:      Number of harmonic components to compute (before
                            mask applied)
        fmin:               Lower bound for tune
        fmax:               Upper bound for tune
        hann:               Turn on Hanning window. Default: :py:obj:`False`

    Returns:
        tunes (ndarray):    Array of tunes qx, qy with shape (len(amp), 2)
    """

    def _gen_part(ring, amp, dim, orbit, ld, nturns):
        invx = np.array(
            [
                [1 / np.sqrt(ld["beta"][0]), 0],
                [ld["alpha"][0] / np.sqrt(ld["beta"][0]), np.sqrt(ld["beta"][0])],
            ]
        )

        invy = np.array(
            [
                [1 / np.sqrt(ld["beta"][1]), 0],
                [ld["alpha"][1] / np.sqrt(ld["beta"][1]), np.sqrt(ld["beta"][1])],
            ]
        )
        part = np.array([orbit] * len(amp)).T + 1.0e-6
        part[dim, :] += amp
        part = internal_lpass(ring, part, nturns=nturns)
        sh = part.shape
        partx = np.reshape(part[0, :], (sh[1], sh[3]))
        partxp = np.reshape(part[1, :], (sh[1], sh[3]))
        party = np.reshape(part[2, :], (sh[1], sh[3]))
        partyp = np.reshape(part[3, :], (sh[1], sh[3]))
        px = np.array([np.matmul(invx, [partx[i], partxp[i]]) for i in range(len(amp))])
        py = np.array([np.matmul(invy, [party[i], partyp[i]]) for i in range(len(amp))])
        return px[:, 0, :] - 1j * px[:, 1, :], py[:, 0, :] - 1j * py[:, 1, :]

    l0, bd, _ = linopt6(ring)
    orbit = l0["closed_orbit"]
    tunes = bd["tune"]

    if amp is not None:
        partx, party = _gen_part(ring, amp, dim, orbit, l0, nturns)
        qx = get_tunes_harmonic(partx, **kwargs)
        qy = get_tunes_harmonic(party, **kwargs)
        tunes = np.vstack((qx, qy)).T

    return np.array(tunes)


def detuning(
    ring: Lattice,
    xm: Optional[float] = 0.3e-4,
    ym: Optional[float] = 0.3e-4,
    npoints: Optional[int] = 3,
    nturns: Optional[int] = 512,
    **kwargs,
):
    """Computes the tunes for a sequence of amplitudes

    This function uses :py:func:`tunes_vs_amp` to compute the tunes for
    the specified amplitudes. Then it fits this data and returns
    the detuning coefficiant dQx/Jx, dQy/Jx, dQx/Jy, dQy/Jy and the
    qx, qy arrays versus x, y arrays

    Parameters:
        ring:       Lattice description
        xm:         Maximum x amplitude
        ym:         Maximum y amplitude
        npoints:    Number of points in each plane
        nturns:     Number of turns for tracking
        method:     ``'laskar'`` or ``'fft'``. Default: ``'laskar'``
        num_harmonics:  Number of harmonic components to compute
                       (before mask applied)
        fmin:       Lower bound for tune
        fmax:       Upper bound for tune
        hann:       Turn on Hanning window. Default: :py:obj:`False`

    Returns:
        q0 (ndarray): qx, qy from horizontal and vertical amplitude scans.
                      Fit results with shape (2, 2). Only the fractional part
                      is returned
        q1 (ndarray): Amplitude detuning coefficients with shape (2, 2)
                      ``[[dQx/dJx, dQy/dJx], [dQx/dJy, dQy/dJy]]``
        x (ndarray): x amplitudes (npoints, )
        q_dx (ndarray): qx, qy tunes as a function of x amplitude (npoints, 2)
                        Only the fractional parts are returned
        y (ndarray): y amplitudes (npoints, )
        q_dy (ndarray): qx, qy tunes as a function of y amplitude (npoints, 2)
                        Only the fractional parts are returned
    """
    lindata0, _, _ = linopt6(ring)
    gamma = (1 + lindata0.alpha * lindata0.alpha) / lindata0.beta

    x = np.linspace(-xm, xm, npoints)
    y = np.linspace(-ym, ym, npoints)
    x2 = x * x
    y2 = y * y

    q_dx = tunes_vs_amp(ring, amp=x, dim=0, nturns=nturns, **kwargs)
    q_dy = tunes_vs_amp(ring, amp=y, dim=2, nturns=nturns, **kwargs)
    q_dx, _ = np.modf(q_dx * ring.periodicity)
    q_dy, _ = np.modf(q_dy * ring.periodicity)

    idx = np.isfinite(q_dx[:, 0]) & np.isfinite(q_dx[:, 1])
    idy = np.isfinite(q_dy[:, 0]) & np.isfinite(q_dy[:, 1])

    fx = np.polyfit(x2[idx], q_dx[idx], 1)
    fy = np.polyfit(y2[idy], q_dy[idy], 1)

    q0 = np.array([[fx[1, 0], fx[1, 1]], [fy[1, 0], fy[1, 1]]])
    q1 = np.array(
        [
            [2 * fx[0, 0] / gamma[0], 2 * fx[0, 1] / gamma[0]],
            [2 * fy[0, 0] / gamma[1], 2 * fy[0, 1] / gamma[1]],
        ]
    )

    return q0, q1, x, q_dx, y, q_dy


def chromaticity(
    ring: Lattice,
    method: Optional[str] = "linopt",
    dpm: Optional[float] = 0.02,
    npoints: Optional[int] = 11,
    order: Optional[int] = 3,
    dp: Optional[float] = 0,
    **kwargs,
):
    r"""Computes the non-linear chromaticities

    This function computes the tunes for the specified momentum offsets.
    Then it fits this data and returns the chromaticity up to the given
    order (npoints>order)

    Parameters:
        ring:       Lattice description
        method:     ``'linopt'`` (dfault) returns the tunes from the
          :py:func:`linopt6` function,

          ``'fft'`` tracks a single particle and computes the tunes with fft,

          ``'laskar'`` tracks a single particle and computes the tunes with
          NAFF.
        dpm:        Maximum momentum deviation
        npoints:    Number of momentum deviations
        order:      Order of the fit
        dp:         Central momentum deviation

    Returns:
        fit (ndarray): Horizontal Vertical polynomial coefficients
                       from ``np.polyfit`` of shape (2, order+1)
        dpa: dp array of shape (npoints,)
        qz: tune array of shape (npoints, 2)
    """
    if order == 0:
        return get_chrom(ring, dp=dp, **kwargs)
    elif order > npoints - 1:
        raise ValueError("order should be smaller than npoints-1")
    else:
        dpa = np.linspace(-dpm, dpm, npoints)
        qz = []
        for dpi in dpa:
            qz.append(
                get_tune(ring, method=method, dp=dp + dpi, remove_dc=True, **kwargs)
            )
        fit = np.polyfit(dpa, qz, order)[::-1]
        fitx = fit[:, 0] * factorial(np.arange(order + 1))
        fity = fit[:, 1] * factorial(np.arange(order + 1))
        return np.array([fitx, fity]), dpa, np.array(qz)


def gen_detuning_elem(ring: Lattice, orbit: Optional[Orbit] = None) -> Element:
    """Generates a detuning element

    Parameters:
        ring:           Lattice description
        orbit:          Avoids looking for initial the closed orbit if it is
          already known ((6,) array).

    Returns:
        detuneElem:     Element reproducing the detuning of the ring with
          amplitude and momentum
    """
    if orbit is None:
        orbit, _ = find_orbit(ring)
    lindata0, _, _ = ring.linopt6(get_chrom=False, orbit=orbit)
    xsi = get_chrom(ring.disable_6d(copy=True))
    r0, r1, x, q_dx, y, q_dy = detuning(
        ring.radiation_off(copy=True), xm=1.0e-4, ym=1.0e-4, npoints=3
    )
    nonlin_elem = Element(
        "NonLinear",
        PassMethod="DeltaQPass",
        Betax=lindata0.beta[0],
        Betay=lindata0.beta[1],
        Alphax=lindata0.alpha[0],
        Alphay=lindata0.alpha[1],
        chromx_arr=np.array([xsi[0]]),
        chromy_arr=np.array([xsi[1]]),
        A1=r1[0][0],
        A2=r1[0][1],
        A3=r1[1][1],
        T1=-orbit,
        T2=orbit,
        chrom_maxorder=1
    )
    return nonlin_elem
