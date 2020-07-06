"""
Function to compute quantities related to non-linear optics
"""
import numpy
from math import sqrt, atan2
from scipy.special import factorial
from at.lattice import Lattice, check_radiation, uint32_refpts, get_s_pos, \
    bool_refpts
from at.tracking import lattice_pass
from at.physics import HarmonicAnalysis, get_tune, linopt, find_orbit4, \
    get_tunes_harmonic, get_chrom
from at.lattice import Element

__all__ = ['detuning', 'chromaticity']


def tunes_vs_amp(ring, amp=None, dim=0,
                 dp=0, window=1, nturns=512):
    """
    Generates a range of tunes for varying x, y, or z amplitudes
    """
    def _gen_part(ring, amp, dim, orbit, nturns):
        part = numpy.array([orbit, ] * len(amp)).T + 1.0e-6
        part[dim, :] += amp
        part = lattice_pass(ring, part, nturns=nturns)
        sh = part.shape
        partx = numpy.reshape(part[0, :], (sh[1], sh[3]))
        partxp = numpy.reshape(part[1, :], (sh[1], sh[3]))
        party = numpy.reshape(part[2, :], (sh[1], sh[3]))
        partyp = numpy.reshape(part[3, :], (sh[1], sh[3]))
        return partx+1j*partxp, party+1j*partyp

    orbit, _ = find_orbit4(ring,dp)
    q0 = get_tune(ring)
    fmin = numpy.maximum([0,0],q0-window)
    fmax = numpy.minimum([1,1],q0+window)
    tunes = []

    if amp is not None:
        partx, party = _gen_part(ring, amp, dim, orbit, nturns)
        qx = get_tunes_harmonic(partx, 'laskar',
                                fmin=fmin[0], fmax=fmax[0])
        qy = get_tunes_harmonic(party, 'laskar',
                                fmin=fmin[1], fmax=fmax[1])
        tunes = numpy.vstack((qx, qy)).T

    return numpy.array(tunes)


def detuning(ring, xm=1.0e-4, ym=1.0e-4, npoints=3, dp=0, window=1, nturns=512):
    """
    This function uses tunes_vs_amp to compute the tunes for
    the specified amplitudes. Then it fits this data and returns
    result for dQx/dx, dQy/dx, dQx/dy, dQy/dy
    """
    lindata0, _, _, _ = linopt(ring, dp=dp)
    gamma = (1 + lindata0.alpha * lindata0.alpha) / lindata0.beta

    x = numpy.linspace(-xm, xm, npoints)
    y = numpy.linspace(-ym, ym, npoints)
    x2 = x * x
    y2 = y * y

    q_dx =  tunes_vs_amp(ring, amp=x, dim=0, dp=dp, window=window, nturns=nturns)
    q_dy =  tunes_vs_amp(ring, amp=y, dim=2, dp=dp, window=window, nturns=nturns)

    fx = numpy.polyfit(x2, q_dx, 1)
    fy = numpy.polyfit(y2, q_dy, 1)

    q0 = [[fx[1, 0], fx[1, 1]], 
          [fy[1, 0], fy[1, 1]]]
    q1 = [[2 * fx[0, 0] / gamma[0],2 * fx[0, 1] / gamma[0]],
          [2 * fy[0, 0] / gamma[1],2 * fy[0, 1] / gamma[1]]]

    return numpy.array(q0), numpy.array(q1), x, q_dx, y, q_dy


def chromaticity(ring, method='linopt', dpm=0.02, npoints=11, order=3, dp=0,**kwargs):
    """
    This function uses tunes_vs_amp to compute the tune for
    the specified amplitudes. Then it fits this data and returns
    the chromaticity up to the given order (npoints>order)
    Q'n = 1/n!*dQ/ddp
    """
    if order ==0:
       return get_chrom(ring,dp=dp,**kwargs)
    elif order>npoints-1:
       raise ValueError('order should be smaller than npoints-1')
    else:    
        dpa = numpy.linspace(-dpm, dpm, npoints)
        qz = []
        for dpi in dpa:
            qz.append(get_tune(ring,method=method,dp=dp+dpi,remove_dc=True,**kwargs))
        fit = numpy.polyfit(dpa, qz, order)[::-1]
        fitx = fit[:,0]/factorial(numpy.arange(order+1))
        fity = fit[:,1]/factorial(numpy.arange(order+1)) 
        return numpy.array([fitx,fity]), dpa, numpy.array(qz)




