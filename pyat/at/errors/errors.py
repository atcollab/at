import numpy
from typing import Callable
import functools
from at.lattice import refpts_iterator, Lattice
from at.lattice import shift_elem, rotate_elem
from at.physics import find_orbit, get_optics
from at.tracking import lattice_pass
from scipy.stats import truncnorm, norm

__all__ = ['find_orbit_err', 'get_optics_err', 'get_ring_with_errors']


_BPM_ATTRS = ['BPMGain', 'BPMOffset', 'BPMTilt']
_ERR_ATTRS = ['PolynomBErr', 'PolynomAErr', 'ShiftErr', 'RotationErr']


def _truncated_randn(shape, sigma=1.0, mean=0.0, truncation=None, seed=None):
    if seed is not None:
        numpy.random.seed(seed)
    npoints = numpy.prod(shape)
    sigmas = numpy.broadcast_to(sigma, shape)
    means = numpy.broadcast_to(mean, shape)
    if truncation is not None:
        pts = truncnorm.rvs(-truncation, truncation, size=npoints)
    else:
        pts = norm.rvs(size=npoints)
    return numpy.reshape(pts, shape) * sigmas + means


def assign_errors(ring, key, truncation=None, seed=None, **kwargs):
    elements = ring.get_elements(key)
    for attr in _BPM_ATTRS+_ERR_ATTRS:
        val = kwargs.pop(attr, None)
        if val is not None:
            vals = numpy.broadcast_to(val, (len(elements),
                                            numpy.size([val])))                                            
            rv = _truncated_randn(vals.shape, sigma=vals,
                                  truncation=truncation, seed=seed)
            [setattr(e, attr, v) for e, v in zip(elements, rv)]
            

def _rotmat(theta):
    cs = numpy.cos(theta)
    sn = numpy.sin(theta)
    return numpy.array([[cs, -sn], [sn, cs]])
    

def _apply_bpm_orbit_error(ring, orbit, refpts):
    for e, o6 in zip(refpts_iterator(ring, refpts), orbit):
        if hasattr(e, 'BPMGain'):
            o6[[0, 2]] *= e.BPMGain
        if hasattr(e, 'BPMOffset'):
            o6[[0, 2]] += e.BPMOffset
        if hasattr(e, 'BPMTilt'):
            o6[[0, 2]] = _rotmat(e.BPMTilt) @ o6[[0, 2]]
    return orbit
    
    
def _apply_bpm_track_error(ring, trajectory, refpts):
    for traj in trajectory.T:
        for e, o6 in zip(refpts_iterator(ring, refpts), traj):
            if hasattr(e, 'BPMGain'):
                o6[:, [0, 2]] *= e.BPMGain
            if hasattr(e, 'BPMOffset'):
                o6[:, [0, 2]] += e.BPMOffset
            if hasattr(e, 'BPMTilt'):
                o6[:, [0, 2]] = [_rotmat(e.BPMTilt) @ o2 
                                 for o2 in o6[:, [0, 2]]]
    return trajectory
    
    
def _apply_alignment_errors(ring):
    refpts = [(hasattr(e, 'ShiftErr') or hasattr(e, 'RotationErr'))
              for e in ring]
    ring = ring.replace(refpts)
    for e in ring[refpts]:
        shift = getattr(e, 'ShiftErr', None) 
        rots = getattr(e, 'RotationsErr', None)   
        if shift is not None:   
            shift_elem(e, shift[0], shift[1])
        if rots is not None:
            rotate_elem(e, tilt=rots[0], pitch=rots[1], yaw=rots[2])
    return ring
    
               
def _apply_field_errors(ring): 

    def sanitize(e):
        pola = getattr(e, 'PolynomA')
        polb = getattr(e, 'PolynomB')
        lm = max(len(pola), len(polb))
        pan = numpy.zeros(lm)
        pbn = numpy.zeros(lm)
        pan[:len(pola)] = pola
        pbn[:len(polb)] = polb
        e.PolynomA = pan
        e.PolynomB = pbn 
        e.MaxOrder = lm-1
 
    def set_polerr(e, pname):
        pol = getattr(e, pname)
        err = getattr(e, pname+'Err')
        pn = numpy.zeros(max(len(pol), len(err)))
        if len(pol) > len(err):
            pn[:len(err)] = numpy.add(pol[:len(err)], err)
            pn[len(err):] = pol[len(err):]
        else:
            pn[:len(pol)] = numpy.add(err[:len(pol)], pol)  
            pn[len(pol):] = err[len(pol):]
        setattr(e, pname, pn)
    
    refpts = [(hasattr(e, 'PolynomBErr') or hasattr(e, 'PolynomAErr'))
              and hasattr(e, 'PolynomB') and hasattr(e, 'PolynomA')
              for e in ring]
    ring = ring.replace(refpts)                            
    for e in ring[refpts]:
        if hasattr(e, 'PolynomBErr'):  
            set_polerr(e, 'PolynomB')
        if hasattr(e, 'PolynomAErr'):        
            set_polerr(e, 'PolynomA') 
        sanitize(e) 
    return ring      


def _orbit_errors(func) -> Callable:
    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs):
        ring = _apply_field_errors(ring)     
        ring = _apply_alignment_errors(ring)     
        orbit0, orbit = func(ring, *args, **kwargs)
        refpts = kwargs.get('refpts', None)
        if refpts is not None:
            orbit = _apply_bpm_orbit_error(ring, orbit, refpts)
        return orbit0, orbit
    return wrapper
    
    
def _linopt_errors(func) -> Callable:
    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs): 
        ring = _apply_field_errors(ring)     
        ring = _apply_alignment_errors(ring)       
        ld0, bd, ld = func(ring, *args, **kwargs)
        refpts = kwargs.get('refpts', None)
        if refpts is not None:
            ld.closed_orbit = \
                _apply_bpm_orbit_error(ring, ld.closed_orbit, refpts)
        return ld0, bd, ld
    return wrapper
    
    
def _track_errors(func) -> Callable:
    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs): 
        ring = _apply_field_errors(ring)     
        ring = _apply_alignment_errors(ring)       
        rout = func(ring, *args, **kwargs)
        refpts = kwargs.get('refpts', None)
        if refpts is not None:
            rout = _apply_bpm_track_error(ring, rout, refpts)
        return rout
    return wrapper    
    

def get_ring_with_errors(ring):
    ring = _apply_field_errors(ring) 
    ring = _apply_alignment_errors(ring) 
    for e in ring:
       for a in _ERR_ATTRS:
           if hasattr(e, a):
               delattr(e, a)
    return ring


find_orbit_err = _orbit_errors(find_orbit)
Lattice.find_orbit_err = find_orbit_err
get_optics_err = _linopt_errors(get_optics)
Lattice.get_optics_err = get_optics_err
lattice_pass_err = _track_errors(lattice_pass)
