import numpy
from typing import Callable
import functools
from at.lattice import bool_refpts, Lattice
from at.physics import find_orbit4, find_sync_orbit
from at.physics import find_orbit6, find_orbit


__all__ = ['find_orbit4_err', 'find_sync_orbit_err',
           'find_orbit6_err', 'find_orbit_err']


def _rotmat(theta):
    cs = numpy.cos(theta)
    sn = numpy.sin(theta)
    return numpy.array([[cs, -sn], [sn, cs]])


def _apply_bpm_errors(func) -> Callable:
    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs):
        orbit0, orbit = func(ring, *args, **kwargs)
        refpts = kwargs.pop('refpts', None)
        if refpts is not None:
            for e, o6 in zip(ring[bool_refpts(refpts, len(ring))], orbit):
                if hasattr(e, 'BPMGain'):
                    o6[[0, 2]] *= e.BPMGain
                if hasattr(e, 'BPMOffset'):
                    o6[[0, 2]] += e.BPMOffset
                if hasattr(e, 'BPMTilt'):
                    o6[[0, 2]] = _rotmat(e.BPMTilt) @ o6[[0, 2]]
        return orbit0, orbit
    return wrapper


find_orbit4_err = _apply_bpm_errors(find_orbit4)
find_sync_orbit_err = _apply_bpm_errors(find_sync_orbit)
find_orbit6_err = _apply_bpm_errors(find_orbit6)
find_orbit_err = _apply_bpm_errors(find_orbit)
Lattice.find_orbit4 = find_orbit4_err
Lattice.find_sync_orbit_err = find_sync_orbit_err
Lattice.find_orbit6_err = find_orbit6_err
Lattice.find_orbit_err = find_orbit_err
