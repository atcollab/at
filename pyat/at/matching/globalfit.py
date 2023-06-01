"""
Collection of functions to fit various global parameters
such as the tune and the chromaticity
"""
from typing import Optional
import numpy
from ..lattice import Lattice, Refpts, set_value_refpts
from ..lattice import AtWarning, AtError
from ..physics import get_tune, get_chrom


__all__ = ['fit_tune', 'fit_chrom']


def _get_tune(ring: Lattice, dp: float, **kwargs):
    get_integer = kwargs.pop('fit_integer', False)
    return get_tune(ring, dp=dp, get_integer=get_integer)[0:2]


def _get_chrom(ring: Lattice, dp: float, **kwargs):
    return get_chrom(ring, dp=dp)[0:2]


def _fit_tune_chrom(ring: Lattice, index: int, func,
                    refpts1: Refpts, refpts2: Refpts, newval,
                    tol: Optional[float] = 1.0e-12,
                    dp: Optional[float] = 0, niter: Optional[int] = 3,
                    regex=False, **kwargs):

    def _get_resp(ring: Lattice, index: int, func, refpts, attname,
                  delta, dp, regex=False, **kwargs):
        set_value_refpts(ring, refpts, attname, delta, index=index,
                         increment=True, regex=regex)
        datap = func(ring, dp, **kwargs)
        set_value_refpts(ring, refpts, attname, -2*delta, index=index,
                         increment=True, regex=regex)
        datan = func(ring, dp, **kwargs)
        set_value_refpts(ring, refpts, attname, delta, index=index,
                         increment=True, regex=regex)
        data = numpy.subtract(datap, datan)/(2*delta)
        return data

    def _fit(ring, index, func, refpts1, refpts2, newval, J,
             dp: Optional[float] = 0, regex=False, **kwargs):
        val = func(ring, dp, **kwargs)
        dk = numpy.linalg.solve(J, numpy.subtract(newval, val))
        set_value_refpts(ring, refpts1, 'PolynomB', dk[0], index=index,
                         increment=True, regex=regex)
        set_value_refpts(ring, refpts2, 'PolynomB', dk[1], index=index,
                         increment=True, regex=regex)
        val = func(ring, dp, **kwargs)
        sumsq = numpy.sum(numpy.square(numpy.subtract(val, newval)))
        return sumsq

    delta = 1e-6 * 10 ** index
    dq1 = _get_resp(ring, index, func, refpts1, 'PolynomB',
                    delta, dp, regex=regex, **kwargs)
    dq2 = _get_resp(ring, index, func, refpts2, 'PolynomB',
                    delta, dp, regex=regex, **kwargs)
    J = [[dq1[0], dq2[0]], [dq1[1], dq2[1]]]

    n = 0
    sumsq = tol+1
    print('Initial value', func(ring, dp, **kwargs))
    while sumsq > tol and n < niter:
        sumsq = _fit(ring, index, func, refpts1, refpts2, newval,
                     J, dp=dp, regex=regex, **kwargs)
        print('iter#', n, 'Res.', sumsq)
        n += 1
    print('Final value', func(ring, dp, **kwargs), '\n')
    return


def fit_tune(ring: Lattice, refpts1: Refpts, refpts2: Refpts, newval,
             tol: float = 1.0e-12,
             dp: Optional[float] = 0, niter: int = 3, regex=False,
             **kwargs) -> None:
    """Fits the tunes using 2 families

    Args:
        ring:       Lattice description
        refpts1:    Selection of the 1st family
        refpts2:    Selection of the 2nd family
        newval:     New tunes, in case an non-zero integer part
                    is provided, fit_integer is set to True

    Keyword arguments:
        tol:        Tolerance for the matching; Default: 1.0e-12
        dp:         Momentum deviation. Default: 0
        niter:      Maximum number of iterations. Default 3
        fit_integer: bool (default=False), use integer tune
        regex:      Using regular expressions for refpt string matching;
                    Default: False

    Typical usage:
    at.fit_tune(ring, refpts1, refpts2, [0.1,0.25])
    """
    print('\nFitting Tune...')
    if numpy.any(numpy.floor(newval) != 0.0):
        kwargs['fit_integer'] = True
    _fit_tune_chrom(ring, 1, _get_tune, refpts1, refpts2, newval, tol=tol,
                    dp=dp, niter=niter, regex=regex, **kwargs)


def fit_chrom(ring: Lattice, refpts1: Refpts, refpts2: Refpts, newval,
              tol: Optional[float] = 1.0e-12,
              dp: Optional[float] = 0, niter: Optional[int] = 3, regex=False,
              **kwargs) -> None:
    """Fit the chromaticities using 2 families

    Args:
        ring:       Lattice description
        refpts1:    Selection of the 1st family
        refpts2:    Selection of the 2nd family
        newval:     New tunes

    Keyword arguments:
        tol:        Tolerance for the matching; Default: 1.0e-12
        dp:         Momentum deviation. Default: 0
        niter:      Maximum number of iterations. Default 3
        regex:      Using regular expressions for refpt string matching;
                    Default: False

    Typical usage:
    at.fit_chrom(ring, refpts1, refpts2, [10,5])
    """
    print('\nFitting Chromaticity...')
    _fit_tune_chrom(ring, 2, _get_chrom, refpts1, refpts2, newval, tol=tol,
                    dp=dp, niter=niter, regex=regex)
