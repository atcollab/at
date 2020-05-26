"""
Collection of functions to fit various global parameters
such as the tune and the chromaticity
"""
import numpy
from at.lattice import get_refpts, set_value_refpts
from at.physics import linopt

__all__ = ['fit_tune', 'fit_chrom']


def _get_tune(ring, dp=0):
    _, tune, _, _ = linopt(ring, dp=dp)
    return tune


def _get_chrom(ring, dp=0):
    _, _, chrom, _ = linopt(ring, dp=dp, get_chrom=True)
    return chrom


def _fit_tune_chrom(ring, order, func, key1, key2, newval, tol=1.0e-12,
                    dp=0, niter=3):

    def _get_resp(ring, order, func, keyi, attname, delta, dp=0):
        set_value_refpts(ring, keyi, attname, delta, order=order,
                         increment=True)
        datap = func(ring, dp=dp)
        set_value_refpts(ring, keyi, attname, -2*delta, order=order,
                         increment=True)
        datan = func(ring, dp=dp)
        set_value_refpts(ring, keyi, attname, delta, order=order,
                         increment=True)
        data = numpy.subtract(datap, datan)/(2*delta)
        return data

    fam1i = get_refpts(ring, key1)
    fam2i = get_refpts(ring, key2)
    try:
        assert fam1i.size != 0 and fam2i.size != 0
    except AssertionError:
        raise ValueError('The selected elements are not found in ring')

    delta = 1e-6*10**(order)
    val = func(ring, dp=dp)
    dq1 = _get_resp(ring, order, func, fam1i, 'PolynomB', delta, dp=dp)
    dq2 = _get_resp(ring, order, func, fam2i, 'PolynomB', delta, dp=dp)
    J = [[dq1[0], dq2[0]], [dq1[1], dq2[1]]]
    dk = numpy.linalg.solve(J, numpy.subtract(newval, val))
    set_value_refpts(ring, fam1i, 'PolynomB', dk[0], order=order,
                     increment=True)
    set_value_refpts(ring, fam2i, 'PolynomB', dk[1], order=order,
                     increment=True)

    val = func(ring, dp=dp)
    sumsq = numpy.sum(numpy.square(numpy.subtract(val, newval)))
    if sumsq > tol and niter > 0:
        _fit_tune_chrom(ring, order, func, key1, key2, newval, tol=tol,
                        dp=dp, niter=niter-1)
    else:
        return


def fit_tune(ring, key1, key2, newval, tol=1.0e-12, dp=0, niter=3):
    """
    Function to fit the tune of the ring, using 2 families famname1
    and famname2 are used to select elements based on wildcards
    compared to FamName

    Args:
        ring: lattice for which the tune or chromaticity needs to be matched
        key1/2: can be:
             1) an element instance, will return all elements of the same type
                in the lattice, e.g. key=Drift('d1', 1.0)
             2) an element type, will return all elements of that type in the
                lattice, e.g. key=at.elements.Sextupole
             3) a string to match against elements' FamName, supports Unix
                shell-style wildcards, e.g. key='BPM_*1'
        newval: new tunes or chromaticities
        tol: tolerance for the matching [default=1.0e-12]
        dp: dp/p at which the values need to be matched [default=0]
        niter: maximum number of iterations to reach tol [default=3]

    Typical usage:
    at.matching.fit_tune(ring, 'QF1*', 'QD2*', [0.1,0.25])
    """
    _fit_tune_chrom(ring, 1, _get_tune, key1, key2, newval, tol=tol, dp=dp,
                    niter=niter)


def fit_chrom(ring, key1, key2, newval, tol=1.0e-12, dp=0, niter=3):
    """
    Function to fit the chromaticity of the ring, using 2 families
    famname1 and famname2 are used to select elements based on wildcards
    compared to FamName

    Args:
        ring: lattice for which the tune or chromaticity needs to be matched
        key1/2: can be:
             1) an element instance, will return all elements of the same type
                in the lattice, e.g. key=Drift('d1', 1.0)
             2) an element type, will return all elements of that type in the
                lattice, e.g. key=at.elements.Sextupole
             3) a string to match against elements' FamName, supports Unix
                shell-style wildcards, e.g. key='BPM_*1'
        newval: new tunes or chromaticities
        tol: tolerance for the matching [default=1.0e-12]
        dp: dp/p at which the values need to be matched [default=0]
        niter: maximum number of iterations to reach tol [default=3]

    Typical usage:
    at.matching.fit_chrom(ring, 'SD*', 'SF*', [10,5])
    """
    _fit_tune_chrom(ring, 2, _get_chrom, key1, key2, newval, tol=tol, dp=dp,
                    niter=niter)
