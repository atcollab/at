"""
Collection of functions to fit various global parameters
such as the tune and the chromaticity
"""
import numpy
from at.lattice import set_value_refpts
from at.physics import linopt

__all__ = ['fit_tune', 'fit_chrom']


def _get_tune(ring, dp=0):
    _, tune, _, _ = linopt(ring, dp=dp)
    return tune


def _get_chrom(ring, dp=0):
    _, _, chrom, _ = linopt(ring, dp=dp, get_chrom=True)
    return chrom


def _fit_tune_chrom(ring, order, func, refpts1, refpts2, newval, tol=1.0e-12,
                    dp=0, niter=3):

    def _get_resp(ring, order, func, refpts, attname, delta, dp=0):
        set_value_refpts(ring, refpts, attname, delta, order=order,
                         increment=True)
        datap = func(ring, dp=dp)
        set_value_refpts(ring, refpts, attname, -2*delta, order=order,
                         increment=True)
        datan = func(ring, dp=dp)
        set_value_refpts(ring, refpts, attname, delta, order=order,
                         increment=True)
        data = numpy.subtract(datap, datan)/(2*delta)
        return data

    delta = 1e-6*10**(order)
    val = func(ring, dp=dp)
    dq1 = _get_resp(ring, order, func, refpts1, 'PolynomB', delta, dp=dp)
    dq2 = _get_resp(ring, order, func, refpts2, 'PolynomB', delta, dp=dp)
    J = [[dq1[0], dq2[0]], [dq1[1], dq2[1]]]
    dk = numpy.linalg.solve(J, numpy.subtract(newval, val))
    set_value_refpts(ring, refpts1, 'PolynomB', dk[0], order=order,
                     increment=True)
    set_value_refpts(ring, refpts2, 'PolynomB', dk[1], order=order,
                     increment=True)

    val = func(ring, dp=dp)
    sumsq = numpy.sum(numpy.square(numpy.subtract(val, newval)))
    if sumsq > tol and niter > 0:
        _fit_tune_chrom(ring, order, func, refpts1, refpts2, newval, tol=tol,
                        dp=dp, niter=niter-1)
    else:
        return


def fit_tune(ring, refpts1, refpts2, newval, tol=1.0e-12, dp=0, niter=3):
    """
    Function to fit the tune of the ring, using 2 families defined by
    refpts

    Args:
        ring: lattice for which the tune needs to be matched
        refpts1/2: refpts for the 2 families
        newval: new tunes
        tol: tolerance for the matching [default=1.0e-12]
        dp: dp/p at which the values need to be matched [default=0]
        niter: maximum number of iterations to reach tol [default=3]

    Typical usage:
    at.matching.fit_tune(ring, refpts1, refpts2, [0.1,0.25])
    """
    _fit_tune_chrom(ring, 1, _get_tune, refpts1, refpts2, newval, tol=tol,
                    dp=dp, niter=niter)


def fit_chrom(ring, refpts1, refpts2, newval, tol=1.0e-12, dp=0, niter=3):
    """
    Function to fit the chromaticity of the ring, using 2 families
    defined by refpts

    Args:
        ring: lattice for which the chromaticity needs to be matched
        ring: lattice for which the chromaticity needs to be matched
        refpts1/2: refpts for the 2 families
        newval: new tunes
        tol: tolerance for the matching [default=1.0e-12]
        dp: dp/p at which the values need to be matched [default=0]
        niter: maximum number of iterations to reach tol [default=3]

    Typical usage:
    at.matching.fit_chrom(ring, refpts1, refpts2, [10,5])
    """
    _fit_tune_chrom(ring, 2, _get_chrom, refpts1, refpts2, newval, tol=tol,
                    dp=dp, niter=niter)
