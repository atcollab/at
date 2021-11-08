"""
Function to compute quantities related to non-linear optics
"""
import numpy
from math import sqrt, atan2
from scipy.special import factorial
from at.lattice import Lattice, check_radiation, uint32_refpts, get_s_pos, \
    bool_refpts
from at.tracking import lattice_pass
from at.physics import HarmonicAnalysis, get_tune, linopt6, find_orbit, \
    get_tunes_harmonic, get_chrom
from at.lattice import Element

__all__ = ['detuning', 'chromaticity', 'gen_detuning_elem']


def tunes_vs_amp(ring, amp=None, dim=0, window=1, nturns=512):
    """
    Generates a range of tunes for varying x, y, or z amplitudes
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
        qx = get_tunes_harmonic(partx, 'laskar')
        qy = get_tunes_harmonic(party, 'laskar')
        tunes = numpy.vstack((qx, qy)).T

    return numpy.array(tunes)


def detuning(ring, xm=0.3e-4, ym=0.3e-4, npoints=3, window=1,
             nturns=512):
    """
    This function uses tunes_vs_amp to compute the tunes for
    the specified amplitudes. Then it fits this data and returns
    result for dQx/dx, dQy/dx, dQx/dy, dQy/dy
    """
    lindata0, _, _ = linopt6(ring)
    gamma = (1 + lindata0.alpha * lindata0.alpha) / lindata0.beta

    x = numpy.linspace(-xm, xm, npoints)
    y = numpy.linspace(-ym, ym, npoints)
    x2 = x * x
    y2 = y * y

    q_dx = tunes_vs_amp(ring, amp=x, dim=0, window=window,
                        nturns=nturns)
    q_dy = tunes_vs_amp(ring, amp=y, dim=2, window=window,
                        nturns=nturns)

    fx = numpy.polyfit(x2, q_dx, 1)
    fy = numpy.polyfit(y2, q_dy, 1)

    q0 = [[fx[1, 0], fx[1, 1]],
          [fy[1, 0], fy[1, 1]]]
    q1 = [[2 * fx[0, 0] / gamma[0], 2 * fx[0, 1] / gamma[0]],
          [2 * fy[0, 0] / gamma[1], 2 * fy[0, 1] / gamma[1]]]

    return numpy.array(q0), numpy.array(q1), x, q_dx, y, q_dy


def chromaticity(ring, method='linopt', dpm=0.02, npoints=11, order=3, dp=0,
                 **kwargs):
    """
    This function uses computes the tunes to compute the tune for
    the specified momentum offsets. Then it fits this data and returns
    the chromaticity up to the given order (npoints>order)
    OUTPUT
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


def gen_detuning_elem(ring, orbit=None):
    """
    Generates an element that for detuning with amplitude
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
