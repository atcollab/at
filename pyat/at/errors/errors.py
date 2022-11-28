import numpy
from typing import Callable
import functools
from at.lattice import refpts_iterator, Lattice
from at.physics import find_orbit4, find_sync_orbit
from at.physics import find_orbit6, find_orbit
from at.physics import linopt2, linopt4, linopt6
from at.physics import get_optics


__all__ = ['find_orbit4_err', 'find_sync_orbit_err',
           'find_orbit6_err', 'find_orbit_err']


def _rotmat(theta):
    cs = numpy.cos(theta)
    sn = numpy.sin(theta)
    return numpy.array([[cs, -sn], [sn, cs]])
    

def _apply_bpm_error(ring, orbit, refpts):
    if refpts is not None:
        for e, o6 in zip(refpts_iterator(ring, refpts), orbit):
            if hasattr(e, 'BPMGain'):
                o6[[0, 2]] *= e.BPMGain
            if hasattr(e, 'BPMOffset'):
                o6[[0, 2]] += e.BPMOffset
            if hasattr(e, 'BPMTilt'):
                o6[[0, 2]] = _rotmat(e.BPMTilt) @ o6[[0, 2]]
    return orbit
 
    
def _apply_field_errors(ring): 

    def sanitize(e):
        pola = getattr(e, 'PolynomA', None)
        polb = getattr(e, 'PolynomB', None)
        lm = max(len(pola), len(polb))
        pan = numpy.zeros(lm)
        pbn = numpy.zeros(lm)
        pan[:len(pola)] = pola
        pbn[:len(polb)] = polb
        e.PolynomA = pan
        e.PolynomB = pbn 
        e.MaxOrder = lm-1
 
    def set_polerr(e, pname):
        pol = getattr(e, pname, None)
        err = getattr(e, pname+'Err', None)
        pn = numpy.zeros(max(len(pol), len(err)))
        if len(pol) > len(err):
            pn[:len(err)] = numpy.add(pol[:len(err)], err)
            pn[len(err):] = pol[len(err):]
        else:
            pn[:len(pol)] = numpy.add(err[:len(pol)], pol)  
            pn[len(pol):] = err[len(pol):]
        setattr(e, pname, pn)
    
    refpts = [(hasattr(e, 'PolynomBErr') 
              or hasattr(e, 'PolynomAErr'))
              and hasattr(e, 'PolynomB') 
              and hasattr(e, 'PolynomA')
              for e in ring]
    ring = ring.replace(refpts)                            
    for e in ring[refpts]:
        if hasattr(e, 'PolynomBErr'):  
            set_polerr(e, 'PolynomB')
        if hasattr(e, 'PolynomAErr'):        
            set_polerr(e, 'PolynomA') 
        sanitize(e) 
    return ring      


def _apply_orbit_errors(func) -> Callable:
    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs):
        ring = _apply_field_errors(ring)          
        orbit0, orbit = func(ring, *args, **kwargs)
        refpts = kwargs.get('refpts', None)
        if refpts is not None:
            orbit = _apply_bpm_error(ring, orbit, refpts)
        return orbit0, orbit
    return wrapper
    
    
def _apply_linopt_errors(func) -> Callable:
    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs): 
        ring = _apply_field_errors(ring)          
        ld0, bd, ld = func(ring, *args, **kwargs)
        refpts = kwargs.get('refpts', None)
        if refpts is not None:
            ld.closed_orbit = _apply_bpm_error(ring, 
                                               ld.closed_orbit,
                                               refpts)
        return ld0, bd, ld
    return wrapper


find_orbit4_err = _apply_orbit_errors(find_orbit4)
find_sync_orbit_err = _apply_orbit_errors(find_sync_orbit)
find_orbit6_err = _apply_orbit_errors(find_orbit6)
find_orbit_err = _apply_orbit_errors(find_orbit)
Lattice.find_orbit4 = find_orbit4_err
Lattice.find_sync_orbit_err = find_sync_orbit_err
Lattice.find_orbit6_err = find_orbit6_err
Lattice.find_orbit_err = find_orbit_err
linopt2_err = _apply_linopt_errors(linopt2)
linopt4_err = _apply_linopt_errors(linopt4)
linopt6_err = _apply_linopt_errors(linopt6)
get_optics_err = _apply_linopt_errors(get_optics)
Lattice.linopt2_err = linopt2_err
Lattice.linopt4_err = linopt4_err
Lattice.linopt6_err = linopt6_err
Lattice.get_optics_err = get_optics_err
