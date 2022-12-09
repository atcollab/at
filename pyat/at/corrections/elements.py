"""
Classes for response matrix observables and variables
"""

from typing import List
import numpy
from at.physics import find_orbit, get_optics
from at.tracking import lattice_pass
from fnmatch import fnmatch


_DATA_NAMES = {'trajectory': ['X', 'XP', 'Y', 'YP', 'DP', 'CT'],
               'closed_orbit': ['X', 'XP', 'Y', 'YP', 'DP', 'CT'],
               'beta': ['BETX', 'BETY'],
               'alpha': ['ALFX', 'ALFY'],
               'mu': ['MUX', 'MUY'],
               'dispersion': ['DX', 'DPX', 'DY', 'DPY'],
               'W': ['WX', 'WY']}


def sum_polab(ring, refpts, attrname, index):
    k = [getattr(ring[r], attrname)[index] *
         getattr(ring[r], 'Length') for r in refpts]
    return numpy.sum(k)


class Variable(object):

    def __init__(self, name, delta, refpts=None, attname=None, index=None,
                 getsetf=None, args=(), kwargs={}):
        self.name = name
        self.attname = attname
        self.refpts = refpts
        self.delta = delta
        self.attindex = index
        self.args = args
        self.kwargs = kwargs
        self.ufun = False
        self.index = index
        if getsetf is None:
            assert refpts is not None, \
                'repfts required for element variable'
            assert attname is not None, \
                'attname required for element variable'
            if index is None:
                self.setf = setattr
                self.getf = getattr
            else:
                self.setf = self._setf
                self.getf = self._getf
        else:
            self.ufun = True
            self.setf, self.getf = getsetf
        super(Variable, self).__init__()

    def _setf(self, elem, attrname, value):
        getattr(elem, attrname)[self.index] = value

    def _getf(self, elem, attrname):
        return getattr(elem, attrname)[self.index]

    def set(self, ring, value):
        if self.ufun:
            self.setfun(ring, value, *self.args, **self.kwargs)
        else:
            for elem in ring.select(self.refpts):
                self.setf(elem, self.attname, value)

    def get(self, ring):
        if self.ufun:
            return self.getfun(ring, *self.args, **self.kwargs)
        else:
            values = numpy.array([self.getf(elem, self.attname)
                                  for elem in ring.select(self.refpts)])
            return numpy.average(values)


class Observable(object):

    def __init__(self, name, refpts=None, index=None, weight=1, fun=None,
                 args=(), kwargs={}):
        self.refpts = refpts
        self.name = name
        self.parname = None
        self.weight = weight
        self.parindex = index
        self.args = args
        self.kwargs = kwargs

        if fun is None:
            assert index is not None, \
                'index required for orbit or linopt observable'
            assert refpts is not None, \
                'repfts required for orbit or linopt observable'
            assert name in _DATA_NAMES.keys(),\
                '{} observable not defined in at.get_optics'.format(name)
            self.parname = name
            self.name = _DATA_NAMES[name][index]
            if self.parname == 'trajectory':
                fun = lattice_pass
            elif self.parname == 'closed_orbit':
                fun = find_orbit
            else:
                fun = get_optics
        self.fun = fun
        self.funname = fun.__name__

        super(Observable, self).__init__()

    def value(self, ring):
        if self.fun is lattice_pass:
            orb0 = self.kwargs.get('init_coord', numpy.zeros((6, )))
            turn = self.kwargs.get('nturns', 1)
            _ = numpy.squeeze(self.fun(ring, orb0, nturns=self.turn,
                                       refpts=self.refpts))
            val = orb0[self.parindex]
        elif self.fun is find_orbit:
            _, orbit = self.fun(ring, refpts=self.refpts, orbit=orb0)
            val = orbit[:, self.parindex][0]
        elif self.fun is get_optics:
            get_w = (self.parname == 'W')
            twiss_in = self.kwargs.get('twiss_in', None)
            _, _, ld = self.fun(ring, refpts=self.refpts,
                                get_w=get_w, twiss_in=twiss_in)
            val = ld[self.parname][self.parindex][0]
        else:
            val = self.fun(ring, *self.args, **self.kwargs)
        assert numpy.isscalar(val), \
            'Observable return value has to be a scalar'
        return val


class RMElements(List):

    def __init__(self, elements=[], **kwargs):
        self.kwargs = kwargs
        super(RMElements, self).__init__(elements)

    def add_elements(self, ring, elements):
        self += numpy.atleast_1d(elements)


class RMVariables(RMElements):

    def __init__(self, variables=[], **kwargs):
        super(RMVariables, self).__init__(elements=variables, **kwargs)

    def add_elements_refpts(self, ring, name, delta, refpts,
                            attname, index=None):
        attname = numpy.broadcast_to(attname, len(refpts))
        name = numpy.broadcast_to(name, len(refpts))
        delta = numpy.broadcast_to(delta, len(refpts))
        index = numpy.broadcast_to(index, len(refpts))
        self += [Variable(n, d, refpts=r, attname=a, index=i, getsetf=None)
                 for r, a, n, d, i in zip(refpts, attname, name, delta, index)]

    def get(self, ring, mask=None):
        if mask is None:
            mask = numpy.ones(len(self), dtype=bool)
        return [v.get(ring) for v in numpy.array(self)[mask]]

    def set(self, ring, values, mask=None):
        if mask is None:
            mask = numpy.ones(len(self), dtype=bool)
        [v.set(ring, val) for v, val in zip(numpy.array(self)[mask], values)]


class RMObservables(RMElements):

    def __init__(self, observables=[], **kwargs):
        super(RMObservables, self).__init__(elements=observables, **kwargs)

    def add_elements_refpts(self, ring, name, refpts, index=None, weight=1):
        traj = name == 'trajectory'
        name = numpy.broadcast_to(name, len(refpts))
        index = numpy.broadcast_to(index, len(refpts))
        weight = numpy.broadcast_to(weight, len(refpts))
        index = numpy.broadcast_to(index, len(refpts))
        if traj:
            turn = self.kwargs.get('nturns', 1)
            for t in range(turn):
                kwargs = self.kwargs.copy()
                kwargs['nturns'] = t + 1
                self += [Observable(n, r, index=i, weight=w, fun=None,
                                    kwargs=kwargs)
                         for r, n, i, w in zip(refpts, name, index, weight)]
        else:
            self += [Observable(n, r, index=i, weight=w, fun=None,
                                kwargs=self.kwargs)
                     for r, n, i, w in zip(refpts, name, index, weight)]

    def values(self, ring, mask=None):
        if mask is None:
            mask = numpy.ones(len(self), dtype=bool)
        allref = [obs.refpts for obs in numpy.array(self)[mask]
                  if obs.refpts is not None]
        allref = numpy.sort(numpy.unique(allref))
        all_fun = numpy.unique([obs.funname for obs in self])
        trajectory = 'lattice_pass' in all_fun
        linopt = 'get_optics' in all_fun
        orbit = 'find_orbit' in all_fun
        assert not (trajectory and orbit), \
            'Trajectory and orbit observables cannot be combined'
        traj = []
        valref = []
        if trajectory:
            turn = self.kwargs.get('nturns', 1)
            o0 = self.kwargs.get('init_coord', numpy.zeros((6, )))
            traj = lattice_pass(ring, o0, nturns=turn, refpts=allref)
        elif linopt:
            get_w = numpy.any([fnmatch(obs.name, 'W[XY]')
                               for obs in numpy.array(self)[mask]])
            twiss_in = self.kwargs.get('twiss_in', None)
            _, _, valref = get_optics(ring, refpts=allref, get_w=get_w,
                                      twiss_in=twiss_in)
        elif orbit:
            _, valref = find_orbit(ring, refpts=allref)

        vals = []
        for o in numpy.array(self)[mask]:
            if o.refpts is None:
                vals.append(o.value(ring))
            elif o.parname == 'trajectory':
                t = o.kwargs.get('nturns', 1)-1
                vals.append(traj[o.parindex, 0, allref == o.refpts, t][0])
            elif linopt:
                vals.append(valref[o.parname][allref == o.refpts,
                                              o.parindex][0])
            elif orbit:
                vals.append(valref[allref == o.refpts, o.parindex][0])
        return numpy.array(vals)
