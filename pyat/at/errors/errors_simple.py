import numpy as np
from typing import Union, Callable
import numpy
from at.lattice import refpts_iterator, Lattice, Element, Refpts
from at.lattice import shift_elem, rotate_elem
from scipy.stats import truncnorm, norm
import functools
import at.tracking


_BPM_ATTRS = {'BPMGain': (2,), 'BPMOffset': (2,), 'BPMTilt': (1,)}
_ERR_ATTRS = {'PolynomBErr': None, 'PolynomAErr': None,
              'ScalingPolynomBErr': None, 'ScalingPolynomAErr': None,
              'ShiftErr': (2,), 'RotationErr': None}
_ALL_ATTRS = dict(**_BPM_ATTRS, **_ERR_ATTRS)


def _truncated_randn(truncation=None, **kwargs):
    if truncation is not None:
        return truncnorm.rvs(-truncation, truncation, **kwargs)
    else:
        return norm.rvs(**kwargs)


def _sysrand(syst, rand,
             truncation: float,
             seed: Union[int, np.random.Generator]):
    rv = _truncated_randn(size=rand.shape,
                          truncation=truncation, random_state=seed)
    return syst + rv*rand


class ErrorGenerator(object):
    def __init__(self, seed=None):
        self.seed = seed
        self.gen = np.random.default_rng(seed)

    def assign_errors(self, ring: Lattice, refpts: Refpts, truncation=None, seed=None, **kwargs):
        elements = ring[refpts]
        nelems = len(elements)
        if self.seed != seed:
            self.gen = np.random.default_rng(seed)
            self.seed = seed
        for attr, sz in _ALL_ATTRS.items():
            val = kwargs.pop(attr, None)
            if val is not None:
                if isinstance(val, tuple):
                    syst, rand = val
                    rand = np.atleast_2d(rand)
                    syst = np.atleast_2d(syst)
                else:
                    rand = np.atleast_2d(val)
                    syst = np.zeros(rand.shape)
                if sz is None:
                    szsyst = syst.shape[1:]
                    szrand = rand.shape[1:]
                else:
                    szsyst = szrand = sz
                syst = np.broadcast_to(syst, ((nelems,) + szsyst))
                rand = np.broadcast_to(rand, ((nelems,) + szrand))
                for el, s, r in zip(elements, syst, rand):
                    err = _sysrand(s, r, truncation, self.gen)
                    print(el.FamName, err, attr)
                    setattr(el, attr, err)
                    print(getattr(el, attr))


def _apply_bpm_orbit_errors(ring: Lattice, refpts: Refpts, orbit):

    def _rotmat(theta):
        cs = np.cos(theta)
        sn = np.sin(theta)
        return np.array([[cs, sn], [-sn, cs]])

    if refpts is None:
        return orbit
    for elem, o6 in zip(refpts_iterator(ring, refpts), orbit):
        o6 = o6.reshape((-1, 6))
        if hasattr(elem, 'BPMOffset'):
            o6[:, [0, 2]] -= elem.BPMOffset
        if hasattr(elem, 'BPMTilt'):
            o6[:, [0, 2]] = o6[:, [0, 2]] @ _rotmat(elem.BPMTilt).T
        if hasattr(elem, 'BPMGain'):
            o6[:, [0, 2]] *= elem.BPMGain
    return orbit


def _apply_bpm_track_errors(ring: Lattice, refpts: Refpts, trajectory):
    if refpts is None:
        return trajectory
    for traj in trajectory.T:
        _apply_bpm_orbit_errors(ring, refpts, traj)
    return trajectory


def _apply_alignment_errors(ring: Lattice):
    refpts = [(hasattr(e, 'ShiftErr') or hasattr(e, 'RotationErr'))
              for e in ring]
    ring = ring.replace(refpts)
    for elem in ring[refpts]:
        shift = getattr(elem, 'ShiftErr', None)
        rots = getattr(elem, 'RotationErr', None)
        if shift is not None:
            shift_elem(elem, shift[0], shift[1], relative=True)
        if rots is not None:
            rotate_elem(elem, tilt=rots[0], pitch=rots[1],
                        yaw=rots[2], relative=True)
    return ring


def _apply_field_errors(ring: Lattice, **kwargs):
    """Apply the defined field errors"""

    def buildlist(el: Element, attrname: str, enabled: bool, scale: float,
                  index: Union[int, None], plist: list):
        """Append the errors attributes to the polynomial list"""
        def vmask(v, idx, pl):
            if idx is not None:
                if idx < len(v):
                    t = np.zeros(v.shape)
                    t[idx] = v[idx]
                    pl.append(t)
            elif len(v) > 0:
                pl.append(v)

        error = getattr(el, attrname)
        if error is not None and enabled:
            vmask(scale*error, index, plist)

    def get_pol(el: Element, pname: str, iname: str, **kwargs) -> list:
        """Build a list of all 5 polynomials for A or B"""
        pstatic = pname+'Err'
        pdynamic = 'Scaling' + pstatic
        index = kwargs.pop(iname, None)
        plist = [getattr(el, pname)]
        enabled = True
        buildlist(el, pstatic, enabled, 1.0, index, plist)
        buildlist(el, pdynamic, enabled, el.strength, index, plist)
        return plist

    def set_pol(el: Element, pname: str, plist: list, sz: int) -> None:
        """Sum up all polynomials and set the PolynomA/B attribute"""
        pn = np.zeros(sz)
        for plnm in plist:
            pn[:len(plnm)] += plnm
        setattr(el, pname, pn)

    refpts = [hasattr(e, 'PolynomBErr')
              or hasattr(e, 'ScalingPolynomBErr')
              or hasattr(e, 'PolynomAErr')
              or hasattr(e, 'ScalingPolynomAErr')
              for e in ring]
    ring = ring.replace(refpts)

    for elem in ring[refpts]:
        alist = get_pol(elem, 'PolynomA', 'IndexA', **kwargs)
        blist = get_pol(elem, 'PolynomB', 'IndexB', **kwargs)
        psize = max(max(len(p) for p in alist), max(len(p) for p in blist))
        set_pol(elem, 'PolynomA', alist, psize)
        set_pol(elem, 'PolynomB', blist, psize)
    return ring


def _apply_track_errors(func) -> Callable:
    @functools.wraps(func)
    def wrapper(ring, *args, **kwargs):
        if ring.errors_enabled:
            ring = _apply_field_errors(ring, **kwargs)
            ring = _apply_alignment_errors(ring)
            refpts = kwargs.get('refpts', None)
            rout = func(ring, *args, **kwargs)
            rout = _apply_bpm_track_errors(ring, rout, refpts)
            return rout
    return wrapper


def get_errors_enabled(ring):
    """Get the error status of a lattice"""
    return getattr(ring, 'errors_enabled', False)


def set_errors_enabled(ring, val):
    """Set the error status of a lattice"""
    return setattr(ring, 'errors_enabled', val)


def get_mean_std_err(ring, key, attr, index=0):
    vals = [numpy.atleast_1d(getattr(e, attr, 0.0))[index]
            for e in ring.get_elements(key)]
    return numpy.mean(vals), numpy.std(vals)


Lattice.errors_enabled = property(get_errors_enabled,
                                  set_errors_enabled(),
                                  doc="Error enabled")
