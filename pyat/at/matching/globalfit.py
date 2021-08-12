"""
Collection of functions to fit various global parameters
such as the tune and the chromaticity
"""
import numpy
from at.lattice import set_value_refpts
from at.physics import get_tune, get_chrom
from at.matching import matching 

__all__ = ['fit_tune', 'fit_chrom']


def _get_tune(ring, dp=0):
    return get_tune(ring, dp=0)[0:2]


def _get_chrom(ring, dp=0):
    return get_chrom(ring, dp=0)[0:2]


def _fit_tune_chrom(ring, index, func, refpts1, refpts2, newval, tol=1.0e-12,
                    dp=0, niter=3):

    def _get_resp(ring, index, func, refpts, attname, delta, dp=0):
        set_value_refpts(ring, refpts, attname, delta, index=index,
                         increment=True)
        datap = func(ring, dp=dp)
        set_value_refpts(ring, refpts, attname, -2*delta, index=index,
                         increment=True)
        datan = func(ring, dp=dp)
        set_value_refpts(ring, refpts, attname, delta, index=index,
                         increment=True)
        data = numpy.subtract(datap, datan)/(2*delta)
        return data

    def _fit(ring, index, func, refpts1, refpts2, newval, J, dp=0):
        val = func(ring, dp=dp)
        dk = numpy.linalg.solve(J, numpy.subtract(newval, val))
        set_value_refpts(ring, refpts1, 'PolynomB', dk[0], index=index,
                         increment=True)
        set_value_refpts(ring, refpts2, 'PolynomB', dk[1], index=index,
                         increment=True)
        val = func(ring, dp=dp)
        sumsq = numpy.sum(numpy.square(numpy.subtract(val, newval)))
        return sumsq

    delta = 1e-6*10**(index)
    dq1 = _get_resp(ring, index, func, refpts1, 'PolynomB', delta, dp=dp)
    dq2 = _get_resp(ring, index, func, refpts2, 'PolynomB', delta, dp=dp)
    J = [[dq1[0], dq2[0]], [dq1[1], dq2[1]]]
    
    n=0
    sumsq = tol+1
    print('Initial value', func(ring, dp=dp))
    while sumsq > tol and n < niter:
        sumsq= _fit(ring, index, func, refpts1, refpts2, newval, J, dp=dp)
        print('iter#',n,'Res.',sumsq)
        n +=1
    print('Final value', func(ring, dp=dp),'\n')
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
    print('\nFitting Tune...')
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
    print('\nFitting Chromaticity...')
    _fit_tune_chrom(ring, 2, _get_chrom, refpts1, refpts2, newval, tol=tol,
                    dp=dp, niter=niter)
