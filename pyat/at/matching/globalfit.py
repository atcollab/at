"""
Collection of functions to fit various global parameters
such as the tune and the chromaticity
"""
import numpy
from at.lattice import get_refpts, set_value_refpts
from at.physics import linopt

__all__ = ['fit_tune_chrom']


def _get_tune(ring,dp=0):
    _,tune,_,_ = linopt(ring,dp=dp,get_chrom=True)
    return tune

def _get_chrom(ring,dp=0):
    _,_,chrom,_ = linopt(ring,dp=dp,get_chrom=True)
    return chrom


def fit_tune_chrom(ring,varname,key1,key2,newval,tol=1.0e-12,dp=0,niter=3):
    '''
    Function to fit the tune or chromaticity of the ring, using 2 families
    famname1 and famname2 are used to select elements based on wildcards compared to 
    FamName

    Args:
        ring: lattice for which the tune or chromaticity needs to be matched
        varname: tune or chrom, for matching the tune of the chromaticity
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
    at.matching.fit_tune_chrom(ring, 'tune', 'QF1*', 'QD2*', [0.1,0.25])
    at.matching.fit_tune_chrom(ring, 'chrom', 'SD*', 'SF*', [10,5]) 
    '''

    def _get_opt_data(ring,varname,**kwargs): 
        dp=kwargs.pop('dp',0)
        if varname == 'tune':
            return _get_tune(ring,dp=dp)
        elif varname == 'chrom':
            return _get_chrom(ring,dp=dp)
        else:
            print('Field not yet defined')
            return


    def _get_resp_tune_chrom(ring,vari,varname,delta,obsname,dp=0):
        if obsname=='tune':
            order=1
        elif obsname=='chrom':
            order=2    
        set_value_refpts(ring,vari,varname,delta,order=order,increment=True)
        datap = _get_opt_data(ring,obsname,dp=dp)
        set_value_refpts(ring,vari,varname,-2*delta,order=order,increment=True)
        datan = _get_opt_data(ring,obsname,dp=dp)           
        set_value_refpts(ring,vari,varname,delta,order=order,increment=True)
        data = numpy.subtract(datap,datan)/(2*delta)
        return data
   

    if varname=='tune':
        order=1
    elif varname=='chrom':
        order=2

    fam1i = get_refpts(ring,key1)
    fam2i = get_refpts(ring,key2)
    try:
       assert fam1i.size !=0 and fam2i.size !=0
    except AssertionError:
       raise ValueError('The selected elements are not found in ring, no fitting done')

    delta = 1e-6*10**(order)
    val= _get_opt_data(ring,varname) 
    dq1 = _get_resp_tune_chrom(ring,fam1i,'PolynomB',delta,varname,dp=dp)
    dq2 = _get_resp_tune_chrom(ring,fam2i,'PolynomB',delta,varname,dp=dp)
    J = [[dq1[0],dq2[0]],[dq1[1],dq2[1]]]
    dk = numpy.linalg.solve(J,numpy.subtract(newval,val));
    set_value_refpts(ring,fam1i,'PolynomB',dk[0],order=order,increment=True)
    set_value_refpts(ring,fam2i,'PolynomB',dk[1],order=order,increment=True)

    val= _get_opt_data(ring,varname) 
    sumsq = numpy.sum(numpy.square(numpy.subtract(val,newval)))
    if sumsq>tol and niter>0:
        fit_tune_chrom(ring,varname,key1,key2,newval,tol=tol,dp=dp,niter=niter-1)
    else:
        return
