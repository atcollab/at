"""
Collection of functions to fit various global parameters
such as the tune and the chromaticity
"""
import numpy
from at.lattice import get_elements, set_value_refpts
from at.physics import linopt

__all__ = ['fit_tune_chrom']


def _get_tune(ring,dp=0):
    _,tune,_,_ = linopt(ring,dp=dp,get_chrom=True)
    return tune

def _get_chrom(ring,dp=0):
    _,_,chrom,_ = linopt(ring,dp=dp,get_chrom=True)
    return chrom


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


def fit_tune_chrom(ring,varname,famname1,famname2,newval,tol=1.0e-12,dp=0,niter=8):
    if varname=='tune':
        order=1
    elif varname=='chrom':
        order=2
    _,fam1i = get_elements(ring,famname1,return_refpts=True)
    _,fam2i = get_elements(ring,famname2,return_refpts=True)
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
        fit_tune_chrom(ring,varname,famname1,famname2,newval,tol=tol,dp=dp,niter=niter-1)
    else:
        return
