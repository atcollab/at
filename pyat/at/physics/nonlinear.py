"""
Coupled or non-coupled 4x4 non-linear motion
"""
import numpy
from math import sqrt, atan2, pi
from at.lattice import Lattice, check_radiation, uint32_refpts, get_s_pos, \
    bool_refpts
from at.tracking import lattice_pass
from at.physics import HarmonicAnalysis, get_tune, linopt, find_orbit4, \
    get_tunes_harmonic
from at.lattice import Element

__all__ = ['tunes_vs_amp', 'detuning', 'gen_nonlin_elem']


def tunes_vs_amp(ring, xamp=None, yamp=None, zamp=None,
                 dp=0, window=0.1, nturns=512):
    """
    Generates a range of tunes for varying x, y, or z amplitudes
    """
    def _gen_part(ring, amp, dim, orbit, nturns):
        part = numpy.array([orbit, ] * len(amp)).T + 1.0e-6
        part[dim, :] += amp
        part = lattice_pass(ring, part, nturns=nturns)
        sh = part.shape
        partx = numpy.reshape(part[0, :], (sh[1], sh[3]))
        party = numpy.reshape(part[2, :], (sh[1], sh[3]))
        return partx, party

    orbit, _ = find_orbit4(ring)
    q0 = get_tune(ring)

    xtunes = []
    ytunes = []
    ztunes = []

    if zamp is not None:
        for z in zamp:
            ztunes.append(get_tune(ring, dp=dp+z))

    if xamp is not None:
        partx, party = _gen_part(ring, xamp, 0, orbit, nturns)
        qx = get_tunes_harmonic(partx, 'laskar',
                                fmin=q0[0]-window, fmax=q0[0]+window)
        qy = get_tunes_harmonic(party, 'laskar',
                                fmin=q0[1]-window, fmax=q0[1]+window)
        xtunes = numpy.vstack((qx, qy)).T

    if yamp is not None:
        partx, party = _gen_part(ring, yamp, 2, orbit, nturns)
        qx = get_tunes_harmonic(partx, 'laskar',
                                fmin=q0[0] - window, fmax=q0[0] + window)
        qy = get_tunes_harmonic(party, 'laskar',
                                fmin=q0[1] - window, fmax=q0[1] + window)
        ytunes = numpy.vstack((qx, qy)).T

    return numpy.array(xtunes), numpy.array(ytunes), numpy.array(ztunes)


def detuning(ring, xm=1.0e-4, ym=1.0e-4, npoints=3, dp=0):
    """
    This function uses tunes_vs_amp to compute the tunes for
    the specified amplitudes. Then it fits this data and returns
    result for dQx/dx, dQy/dx, dQx/dy, qQy/dy
    """
    lindata0, _, _, _ = linopt(ring, dp=dp)
    gamma = (1 + lindata0.alpha * lindata0.alpha) / lindata0.beta

    x = numpy.linspace(-xm, xm, npoints)
    y = numpy.linspace(-ym, ym, npoints)
    x2 = x * x
    y2 = y * y

    q_dx, q_dy, _ = tunes_vs_amp(ring, xamp=x, yamp=y)
    fx = numpy.polyfit(x2, q_dx, 1)
    fy = numpy.polyfit(y2, q_dy, 1)

    q0 = [fx[1, 0], fx[1, 1], fy[1, 0], fy[1, 1]]
    q1 = [2 * fx[0, 0] / gamma[0],
          2 * fx[0, 1] / gamma[0],
          2 * fy[0, 0] / gamma[1],
          2 * fy[0, 1] / gamma[1]]

    return numpy.array(q0), numpy.array(q1)


def gen_nonlin_elem(ring):
    """
    Generates an element that for detuning with amplitude
    """
    [lindata0, tunes, xsi, lindata] = ring.linopt(dp=0,
                                                  refpts=len(ring),
                                                  get_chrom=True)

    r0, r1 = detuning(ring, xm=1.0e-4, ym=1.0e-4, npoints=3, dp=0)
    orbit, _ = find_orbit4(ring)

    nonlin_elem = Element('NonLinear', PassMethod='DeltaQPass',
                          Betax=lindata.beta[0][0], Betay=lindata.beta[0][1],
                          Alphax=lindata.alpha[0][0],
                          Alphay=lindata.alpha[0][1],
                          Qpx=xsi[0], Qpy=xsi[1], A1=r1[0], A2=r1[1],
                          A3=r1[3], T1=-orbit, T2=orbit)

    return nonlin_elem


Lattice.tunes_vs_amp = tunes_vs_amp
Lattice.detuning = detuning
Lattice.gen_nonlin_elem = gen_nonlin_elem
