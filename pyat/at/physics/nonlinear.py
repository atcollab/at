"""
Non-linear optics
"""
import numpy
from typing import Optional, Sequence
from scipy.special import factorial
from ..lattice import Element, Lattice, check_collective
from ..tracking import lattice_pass
from .orbit import Orbit, find_orbit
from .linear import get_tune, get_chrom, linopt6
from .harmonic_analysis import get_tunes_harmonic

__all__ = ['detuning', 'chromaticity', 'gen_detuning_elem']


@check_collective(False)
def tunes_vs_amp(ring: Lattice, amp: Optional[Sequence[float]] = None,
                 dim: Optional[int] = 0, window=1,
                 nturns: Optional[int] = 512) -> numpy.ndarray:
    r"""Generates a range of tunes for varying x, y, or z amplitudes

    Parameters:
        ring:               Lattice description
        amp:                Oscillation amplitude. Default: tunes from
          :py:func:`.linopt6`
        dim:                Coordinate index of the oscillation
        nturns:             Number of turns

    Returns:
        tunes (ndarray):    Array of tunes
    """

    def _gen_part(ring, amp, dim, orbit, ld, nturns):
        invx = numpy.array([[1/numpy.sqrt(ld['beta'][0]), 0],
                            [ld['alpha'][0]/numpy.sqrt(ld['beta'][0]),
                            numpy.sqrt(ld['beta'][0])]])

        invy = numpy.array([[1/numpy.sqrt(ld['beta'][1]), 0],
                            [ld['alpha'][1]/numpy.sqrt(ld['beta'][1]),
                            numpy.sqrt(ld['beta'][1])]])
        part = numpy.array([orbit, ] * len(amp)).T + 1.0e-6
        part[dim, :] += amp
        part = lattice_pass(ring, part, nturns=nturns)
        sh = part.shape
        partx = numpy.reshape(part[0, :], (sh[1], sh[3]))
        partxp = numpy.reshape(part[1, :], (sh[1], sh[3]))
        party = numpy.reshape(part[2, :], (sh[1], sh[3]))
        partyp = numpy.reshape(part[3, :], (sh[1], sh[3]))
        px = numpy.array([numpy.matmul(invx, [partx[i], partxp[i]])
                          for i in range(len(amp))])
        py = numpy.array([numpy.matmul(invy, [party[i], partyp[i]])
                          for i in range(len(amp))])
        return px[:, 0, :] - 1j*px[:, 1, :], py[:, 0, :] - 1j*py[:, 1, :]

    l0, bd, _ = linopt6(ring)
    orbit = l0['closed_orbit']
    tunes = bd['tune']

    if amp is not None:
        partx, party = _gen_part(ring, amp, dim, orbit, l0, nturns)
        qx = get_tunes_harmonic(partx, method='laskar')
        qy = get_tunes_harmonic(party, method='laskar')
        tunes = numpy.vstack((qx, qy)).T

    return numpy.array(tunes)


def detuning(ring: Lattice,
             xm: Optional[float] = 0.3e-4, ym: Optional[float] = 0.3e-4,
             npoints: Optional[int] = 3, window: Optional[int] = 1,
             nturns: Optional[int] = 512):
    """Computes the tunes for a sequence of amplitudes

    This function uses :py:func:`tunes_vs_amp` to compute the tunes for
    the specified amplitudes. Then it fits this data and returns
    result for dQx/dx, dQy/dx, dQx/dy, dQy/dy

    Parameters:
        ring:       Lattice description
        xm:         Maximum x amplitude
        ym:         Maximum y amplitude
        npoints:    Number of points in each plane
        window:
        nturns:     Number of turns for tracking

    Returns:
        q0 (ndarray):
        q1 (ndarray):
        x (ndarray):          x amplitudes
        q_dx (ndarray):
        y (ndarray):          y amplitudes
        q_dy (ndarray):
    """
    lindata0, _, _ = linopt6(ring)
    gamma = (1 + lindata0.alpha * lindata0.alpha) / lindata0.beta

    x = numpy.linspace(-xm, xm, npoints)
    y = numpy.linspace(-ym, ym, npoints)
    x2 = x * x
    y2 = y * y

    q_dx = tunes_vs_amp(ring, amp=x, dim=0, window=window, nturns=nturns)
    q_dy = tunes_vs_amp(ring, amp=y, dim=2, window=window, nturns=nturns)

    fx = numpy.polyfit(x2, q_dx, 1)
    fy = numpy.polyfit(y2, q_dy, 1)

    q0 = [[fx[1, 0], fx[1, 1]],
          [fy[1, 0], fy[1, 1]]]
    q1 = [[2 * fx[0, 0] / gamma[0], 2 * fx[0, 1] / gamma[0]],
          [2 * fy[0, 0] / gamma[1], 2 * fy[0, 1] / gamma[1]]]

    return numpy.array(q0), numpy.array(q1), x, q_dx, y, q_dy


def chromaticity(ring: Lattice, method: Optional[str] = 'linopt',
                 dpm: Optional[float] = 0.02, npoints: Optional[int] = 11,
                 order: Optional[int] = 3,
                 dp: Optional[float] = 0,
                 **kwargs):
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
        (fitx, fity), dpa, qz
    """
    if order == 0:
        return get_chrom(ring, dp=dp, **kwargs)
    elif order > npoints - 1:
        raise ValueError('order should be smaller than npoints-1')
    else:
        dpa = numpy.linspace(-dpm, dpm, npoints)
        qz = []
        for dpi in dpa:
            qz.append(get_tune(ring, method=method, dp=dp+dpi, remove_dc=True,
                      **kwargs))
        fit = numpy.polyfit(dpa, qz, order)[::-1]
        fitx = fit[:, 0]/factorial(numpy.arange(order + 1))
        fity = fit[:, 1]/factorial(numpy.arange(order + 1))
        return numpy.array([fitx, fity]), dpa, numpy.array(qz)


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
    xsi = get_chrom(ring.radiation_off(copy=True))
    r0, r1, x, q_dx, y, q_dy = detuning(ring.radiation_off(copy=True),
                                        xm=1.0e-4, ym=1.0e-4, npoints=3)
    nonlin_elem = Element('NonLinear', PassMethod='DeltaQPass',
                          Betax=lindata0.beta[0], Betay=lindata0.beta[1],
                          Alphax=lindata0.alpha[0], Alphay=lindata0.alpha[1],
                          Qpx=xsi[0], Qpy=xsi[1],
                          A1=r1[0][0], A2=r1[0][1], A3=r1[1][1],
                          T1=-orbit, T2=orbit)
    return nonlin_elem
