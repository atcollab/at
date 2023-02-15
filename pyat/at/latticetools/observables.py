import sys
from typing import Optional, Union
if sys.version_info.minor < 9:
    from typing import Callable, Container, Iterable, Tuple
    from typing import AbstractSet as Set
else:
    from collections.abc import Callable, Container, Iterable, Set
    Tuple = tuple
from math import pi
from enum import Enum
from itertools import repeat
from ..lattice import Lattice, Orbit, Refpts, All, End
import numpy as np

RefIndex = Union[int, Tuple[int, ...], slice]


def _selector(select: Container[str]):
    if select is None:
        return lambda obs: True
    else:
        return lambda obs: obs.name in select


class Need(Enum):
    """Defines the computation requirements for an :py:class:`Observable`.
    """
    #:  Specify :py:func:`.find_orbit` computation and provide its *orbit*
    #:  t to the evaluation function
    ORBIT = 1
    #:  Specify :py:func:`.find_m44` or :py:func:`.find_m44` computation and
    #:  provide its *m44* or *m66* output to the evaluation function
    MATRIX = 2
    #:  Specify :py:func:`.get_optics` computation and provide its *ringdata*
    #:  t to the evaluation function
    GLOBALOPTICS = 3
    #:  Specify :py:func:`.get_optics` computation and provide its *elemdata*
    #:  t to the evaluation function
    LOCALOPTICS = 4
    #:  Specify :py:func:`.lattice_pass` computation and provide its *r_out*
    #:  t to the evaluation function
    TRAJECTORY = 5
    #:  Specify :py:func:`.envelope_parameters` computation and provide its
    #:  *params* t to the evaluation function
    EMITTANCE = 6
    #:  Associated with LOCALOPTICS, require local optics computation at all
    #:  points: slower but avoids jumps in phase advance
    ALL_POINTS = 7
    #:  Associated with LOCALOPTICS, require the *get_chrom* keyword
    CHROMATICITY = 8


class Observable(object):
    """Base class for Observables. Can be used for user-defined observables"""
    def __init__(self, fun: Callable, *args,
                 name: Optional[str] = None,
                 target=None,
                 weight=1.0,
                 bounds=(0.0, 0.0),
                 needs: Optional[Set[Need]] = None,
                 **kwargs):
        r"""
        Args:
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            fun:            :ref:`evaluation function <base_eval>`
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]
            needs:          Set of requirements. This selects the data provided
              to the evaluation function. *needs* items are members of the
              :py:class:`Need` enumeration
            *args:          Positional arguments of the evaluation function
            **kwargs:       Keyword arguments for the evaluation function

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        .. _base_eval:
        .. rubric:: User-defined evaluation function

        The general for is:

        :pycode:`value = fun(ring, *data, *args, **kwargs)`

        *data* depends on the *needs* argument, and by default is empty.

        *value* is the value of the observable,

        For user-defined evaluation functions using linear optics data or
        emittance data, it is recommended to use
        :py:class:`LocalOpticsObservable`, :py:obj:`GlobalOpticsObservable`
        or :py:class:`EmittanceObservable`.

        Examples:

            >>> def circumference(ring):
            ...     return ring.get_s_pos(len(ring))[0]
            >>> obs = Observable(circumference)

            Defines an Observable for the ring circumference.

            >>> def momentum_compaction(ring):
            ...     return ring.get_mcf()
            >>> obs = Observable(momentum_compaction)

            Defines an Observable for the momentum compaction factor.
        """
        name = fun.__name__ if name is None else name
        self.fun = fun
        self.needs = needs or set()
        self.name = name
        self.target = target
        self.w = weight
        self.lbound, self.ubound = bounds
        self.args = args
        self.kwargs = kwargs
        self.initial = None
        self._value = None

    def __str__(self):
        return "\n".join((self._header(), self._all_lines()))

    @staticmethod
    def _header():
        """Header line"""
        fstring = "{:<12} {:>16}  {:>16}  {:>16}  {:>16}  {:>16} "
        return fstring.format("location", "Initial", "Actual", "Low bound",
                              "High bound", "residual")

    @staticmethod
    def _line(loc, *items):
        def pitem(v):
            if v is None:
                return " {:>16} ".format("-  ")
            elif isinstance(v, np.ndarray) and v.ndim > 0:
                if v.size == 0:
                    return " {:16} ".format("[]")
                else:
                    return " [{:< 10.4} ...] ".format(v.flat[0])
            else:
                return " {: 16.6} ".format(v)

        its = [pitem(v) for v in items]
        return "    {:<12}".format(loc) + "".join(its)

    def _all_lines(self):
        if self._value is None:
            vnow = None
            vmin = None
            vmax = None
        else:
            vnow = np.asarray(self._value)
            if self.target is None:
                vmin = None
                vmax = None
            else:
                target = np.broadcast_to(self.target, vnow.shape)
                vmin = target + self.lbound
                vmax = target + self.ubound
        values = self._line("", self.initial, vnow, vmin, vmax, self.residual)
        return "\n".join((self.name, values))

    def _setup(self, ring: Lattice):
        """Setup function called wen the observable is added to a list"""
        pass

    def evaluate(self, ring: Lattice, *data, initial: bool = False):
        """Compute and store the value of the observable

        Args:
            ring:       Lattice description
            *data:      Raw data, sent to the evaluation function
            initial:    It :py:obj:`None`, store the result as the initia;
              value
        """
        val = self.fun(ring, *data, *self.args, **self.kwargs)
        if initial:
            self.initial = val
        self._value = val

    def clear(self):
        """Clear the initial value"""
        self.initial = None

    @property
    def value(self):
        """Value of the observable"""
        return self._value

    @property
    def weighted_value(self):
        return self._value / self.w

    @property
    def residual(self):
        """residual, computed as
        :pycode:`residual = ((value-target)/weight)**2`"""
        if self._value is None:
            res = None
        else:
            vnow = np.asarray(self.value)
            if self.target is None:
                res = np.broadcast_to(0.0, vnow.shape)
            else:
                diff = vnow - np.broadcast_to(self.target, vnow.shape)
                lb = diff - self.lbound
                ub = diff - self.ubound
                lb[lb >= 0] = 0
                ub[ub <= 0] = 0
                res = ((lb + ub) / self.w) ** 2
        return res

    @staticmethod
    def _idx(index: RefIndex):
        return slice(None) if (index is None) else index

    @classmethod
    def _arrayaccess(cls, index: RefIndex):
        """Access to array elements"""
        idx = cls._idx(index)

        # noinspection PyUnusedLocal
        def array_element(ring, data):
            return data[idx]
        return array_element

    @classmethod
    def _recordaccess(cls, fieldname: str, index: RefIndex):
        """Access to record elements"""
        idx = cls._idx(index)

        # noinspection PyUnusedLocal
        def fun(ring, data):
            return getattr(data, fieldname)[idx]
        return fun

    @staticmethod
    def _set_name(name, param, index):
        """Compute default observable names"""
        if name is None:
            if index is None or index is Ellipsis:
                subscript = ""
            else:
                subscript = "[{}]".format(index)
            if callable(param):
                base = param.__name__
            else:
                base = param
            return base + subscript
        else:
            return name


class _ElementObservable(Observable):
    """Base class for Observables linked to a position in the lattice"""
    def __init__(self, fun: Callable, refpts: Refpts, *args,
                 name: Optional[str] = None,
                 summary: bool = False,
                 statfun: Optional[Callable] = None, **kwargs):
        r"""
        Args:
            fun:            Evaluation function.
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`
            *args:          Optional positional arguments for a user-defined
              evaluation function
            **kwargs:       Optional keyword arguments for a user-defined
              evaluation function

        Keyword Args:
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.
        """
        name = fun.__name__ if name is None else name
        if statfun:
            summary = True
            name = '{}({})'.format(statfun.__name__, name)

            def modfun(ring, *a):
                return statfun(fun(ring, *a), axis=0)
        else:
            modfun = fun
        super().__init__(modfun, *args, name=name, **kwargs)
        self.summary = summary
        self.refpts = refpts
        self._boolrefs = None
        self._excluded = None
        self._locations = []

    def _all_lines(self):
        if self.summary or self._value is None:
            return super()._all_lines()
        else:
            vnow = self._value
            if self.target is None:
                vmin = repeat(None)
                vmax = repeat(None)
            else:
                target = np.broadcast_to(self.target, vnow.shape)
                vmin = target + self.lbound
                vmax = target + self.ubound
            vini = self.initial
            if vini is None:
                vini = repeat(None)
            viter = zip(self._locations, vini, vnow, vmin, vmax, self.residual)
            values = "\n".join(self._line(*vv) for vv in viter)
            return "\n".join((self.name, values))

    def _setup(self, ring: Lattice):
        boolrefs = ring.get_bool_index(self.refpts)
        excluded = ring.get_bool_index(self._excluded)
        boolrefs &= ~excluded
        self._boolrefs = boolrefs
        locs = [el.FamName for el in ring.select(self._boolrefs[:-1])]
        if boolrefs[-1]:
            locs.append("End")
        self._locations = locs

    @staticmethod
    def _idx(index: RefIndex):
        if isinstance(index, tuple):
            return (slice(None),) + tuple(Observable._idx(i) for i in index)
        else:
            return slice(None), Observable._idx(index)


class OrbitObservable(_ElementObservable):
    """Observes the transfer matrix at selected locations"""
    def __init__(self, refpts: Refpts,
                 axis: Optional[RefIndex] = None,
                 name: Optional[str] = None, **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            axis:           Index in the orbit vector, If :py:obj:`None`,
              the whole vector is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated

        Keyword Args:
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        Example:

            >>> obs = OrbitObservable(Monitor, axis=0)

            Observe the horizontal closed orbit at monitor locations
        """
        name = self._set_name(name, 'orbit', axis)
        fun = self._arrayaccess(axis)
        needs = {Need.ORBIT}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class MatrixObservable(_ElementObservable):
    """Observes the closed orbit at selected locations"""
    def __init__(self, refpts: Refpts,
                 axis: Optional[RefIndex] = Ellipsis,
                 name: Optional[str] = None, **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            axis:           Index in the transfer matrix, If :py:obj:`Ellipsis`,
              the whole matrix is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated

        Keyword Args:
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        Example:

            >>> obs = MatrixObservable(Monitor, axis=(0, 1))

            Observe the transfer matrix from origin to monitor locations and
            extract T[0,1]
        """
        name = self._set_name(name, 'matrix', axis)
        fun = self._arrayaccess(axis)
        needs = {Need.MATRIX}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class _GlobalOpticsObservable(Observable):
    def __init__(self, param: str, *args,
                 index: Optional[RefIndex] = None,
                 name: Optional[str] = None,
                 **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            param:          Optics parameter name (see :py:func:`.get_optics`)
              or user-defined evaluation function called as:
              :pycode:`value = fun(ring, ringdata, *args, **kwargs)` and
              returning the value of the Observable
            index:          Index in the parameter array, If :py:obj:`None`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            *args:          Optional positional arguments for a user-defined
              evaluation function
            **kwargs:       Optional keyword arguments for a user-defined
              evaluation function

        Keyword Args:
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.
        """
        needs = {Need.GLOBALOPTICS}
        name = self._set_name(name, param, index)
        if callable(param):
            fun = param
            needs.add(Need.CHROMATICITY)
        else:
            fun = self._recordaccess(param, index)
            if param == 'chromaticity':
                needs.add(Need.CHROMATICITY)
        super().__init__(fun, *args, needs=needs, name=name, **kwargs)


class LocalOpticsObservable(_ElementObservable):
    """Observe a local optics parameter at selected locations"""
    def __init__(self, refpts: Refpts, param: Union[str, Callable], *args,
                 index: Optional[RefIndex] = Ellipsis,
                 name: Optional[str] = None,
                 use_integer: bool = False,
                 **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            param:          Optics parameter name (see :py:func:`.get_optics`)
              or :ref:`user-defined evaluation function <localoptics_eval>`
            index:          Index in the parameter array, If :py:obj:`Ellipsis`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            use_integer:    For  the *'mu'* parameter, compute the
              phase advance at all points to avoid discontinuities (slower)
            *args:          Optional positional arguments for a user-defined
              evaluation function
            **kwargs:       Optional keyword arguments for a user-defined
              evaluation function

        Keyword Args:
            summary:        Set to :py:obj:`True` if the user-defined
             evaluation function returns a single item (see below)
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        .. _localoptics_eval:
        .. rubric:: User-defined evaluation function

        It is called as:

        :pycode:`value = fun(ring, elemdata, *args, **kwargs)`

        *elemdata* if the output of :py:func:`.get_optics`.

        *value* is the value of the Observable and must have one line per
        refpoint. Alternatively, it may be a single line, but then the
        *summary* keyword must be set to :py:obj:`True`

        Examples:

            >>> obs = LocalOpticsObservable(Monitor, 'beta')

            Observe the beta in both planes at all :py:class:`.Monitor`
            locations

            >>> obs = LocalOpticsObservable(Quadrupole, 'beta', index=1,
            ...                        statfun=np.max)

            Observe the maximum vertical beta in Quadrupoles

            >>> def phase_advance(ring, elemdata):
            ...     mu = elemdata.mu
            ...     return mu[-1] - mu[0]
            >>>
            >>> allobs.append(LocalOpticsObservable([33, 101], phase_advance,
            ...               use_integer=True, summary=True))

            The user-defined evaluation function computes the phase-advance
            between the 1st and last given reference points, here the elements
            33 and 101 of the lattice
        """
        needs = {Need.LOCALOPTICS}
        name = self._set_name(name, param, index)
        if callable(param):
            fun = param
            needs.add(Need.CHROMATICITY)
        else:
            fun = self._recordaccess(param, index)
        if use_integer:
            needs.add(Need.ALL_POINTS)

        super().__init__(fun, refpts, *args, needs=needs, name=name, **kwargs)


class RingObservable(_ElementObservable):
    """Observe an attribute of selected lattice elements"""
    def __init__(self, refpts: Refpts, attrname: str,
                 index: Optional[RefIndex] = Ellipsis,
                 name: Optional[str] = None, **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            refpts:         Elements to be observed
              See ":ref:`Selecting elements in a lattice <refpts>`"
            attrname:       Attribute name
            index:          Index in the attribute array, If :py:obj:`Ellipsis`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated

        Keyword Args:
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`

        Example:

            >>> obs = RingObservable(at.Sextupole, 'KickAngle', index=0,
            ...                      statfun=np.sum)

            Observe the sum of horizontal kicks in Sextupoles
        """

        def fun(ring):
            vals = [get_val(ring, el) for el in ring.select(self._boolrefs)]
            return np.array(vals)

        get_val = Observable._recordaccess(attrname, index)
        name = self._set_name(name, attrname, index)
        super().__init__(fun, refpts, name=name, **kwargs)


class TrajectoryObservable(_ElementObservable):
    """Observe trajectory coordinates at selected locations"""
    def __init__(self, refpts: Refpts,
                 index: Optional[RefIndex] = Ellipsis,
                 name: Optional[str] = None, **kwargs):
        r"""
        Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            index:          Index in the orbit array, If :py:obj:`Ellipsis`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated

        Keyword Args:
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.
        """
        name = self._set_name(name, 'trajectory', index)
        fun = self._arrayaccess(index)
        needs = {Need.TRAJECTORY}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class EmittanceObservable(Observable):
    """Observe emittance-related parameters"""
    def __init__(self, param: str,
                 index: Optional[RefIndex] = None,
                 name: Optional[str] = None,
                 **kwargs):
        r"""
        Args:
            param:          Parameter name (see
              :py:func:`.envelope_parameters`)
            index:          Index in the parameter array, If :py:obj:`None`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated

        Keyword Args:
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]

        Example:

            >>> EmittanceObservable('emittances', index=0)

            Observe the horizontal emittance
        """
        name = self._set_name(name, param, index)
        fun = self._recordaccess(param, index)
        needs = {Need.EMITTANCE}
        super().__init__(fun, needs=needs, name=name, **kwargs)


class ObservableList(list):
    """Handles a list of Observables to be evaluated together"""
    def __init__(self, ring: Lattice, obsiter: Iterable[Observable] = ()):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            ring:       Lattice description
            obsiter:    Iterable of :py:class:`Observable`\ s

        Example:

            >>> obslist = ObservableList(ring)

            Create an empty Observable list

            >>> obslist.append(OrbitObservable(Monitor, index=0))
            >>> obslist.append(GlobalOpticsObservable('tune'))
            >>> obslist.append(EmittanceObservable('emittances', index=0))

            Add observables for horizontal closed orbit at Monitor locations,
            tunes and horizontal emittance

            >>> obslist.evaluate(ring, initial=True)

            Compute the values and mark them as the initial value

            >>> obslist.values
            [array([-3.02189464e-09,  4.50695130e-07,  4.08205818e-07,
                     2.37899969e-08, -1.31783561e-08,  2.47230794e-08,
                     2.95310770e-08, -4.05598110e-07, -4.47398092e-07,
                     2.24850671e-09]),
            array([3.81563019e-01, 8.54376397e-01, 1.09060730e-04]),
            1.320391045951568e-10]

            >>> obslist.get_flat_values({'tune', 'emittances[0]'})
            array([3.815630e-01, 8.543764e-01, 1.090607e-04, 1.320391e-10])

            Get a flattened array of tunes and horizontal emittance
        """
        noref = ring.get_bool_index(None)
        self.orbitrefs = noref
        self.opticsrefs = noref
        self.passrefs = noref
        self.matrixrefs = noref
        self.ring = ring
        self.needs = set()
        super().__init__(self._scan(obs) for obs in obsiter)

    # noinspection PyProtectedMember
    def _update_reflists(self, obs):
        needs = obs.needs
        if Need.ORBIT in needs:
            self.orbitrefs |= obs._boolrefs
        if Need.MATRIX in needs:
            self.matrixrefs |= obs._boolrefs
        if Need.LOCALOPTICS in needs:
            if Need.ALL_POINTS in needs:
                self.opticsrefs = self.ring.get_bool_index(All)
            else:
                self.opticsrefs |= obs._boolrefs
        if Need.TRAJECTORY in needs:
            self.passrefs |= obs._boolrefs

    def _scan(self, obs):
        if not isinstance(obs, Observable):
            raise TypeError("{} is not an Observable".format(obs))
        # noinspection PyProtectedMember
        obs._setup(self.ring)
        self.needs |= obs.needs
        self._update_reflists(obs)
        return obs

    def __iadd__(self, other: "ObservableList"):
        if not isinstance(other, ObservableList):
            mess = "Cannot add a {} to an Observable"
            raise TypeError(mess.format(type(other)))
        if other.ring is not self.ring:
            raise TypeError("Observables must be based on the same lattice")
        self.extend(other)
        return self

    def __add__(self, other):
        nobs = ObservableList(self.ring, self)
        nobs += other
        return nobs

    def append(self, obs: Observable):
        super().append(self._scan(obs))

    def extend(self, obsiter: Iterable[Observable]):
        super().extend(self._scan(obs) for obs in obsiter)

    def insert(self, index: int, obs: Observable):
        super().insert(index, self._scan(obs))

    # noinspection PyProtectedMember
    def __str__(self):
        values = "\n".join(obs._all_lines() for obs in self)
        return "\n".join((Observable._header(), values))

    def evaluate(self, ring: Lattice,
                 r_in: Orbit = None,
                 initial: bool = False):
        r"""Compute all the :py:class:`Observable` values

        Args:
            ring:       Lattice description
            r_in:       Optional coordinate input for trajectory observables
            initial:    If :py:obj:`True`, store the values as *initial values*
        """
        def eval(obs):
            obsneeds = obs.needs
            obsrefs = getattr(obs, '_boolrefs', None)
            data = []
            if Need.ORBIT in obsneeds:
                data.append(orbits[obsrefs[self.orbitrefs]])
            if Need.MATRIX in obsneeds:
                data.append(mxdata[obsrefs[self.matrixrefs]])
            if Need.GLOBALOPTICS in obsneeds:
                data.append(rgdata)
            if Need.LOCALOPTICS in obsneeds:
                data.append(eldata[obsrefs[self.opticsrefs]])
            if Need.TRAJECTORY in obsneeds:
                data.append(trajs[obsrefs[self.passrefs]])
            if Need.EMITTANCE in obsneeds:
                data.append(emdata)
            obs.evaluate(ring, *data, initial=initial)

        trajs = orbits = rgdata = eldata = emdata = mxdata = None
        needs = self.needs

        if Need.TRAJECTORY in needs:
            if r_in is None:
                r_in = np.zeros(6)
            r_out = ring.lattice_pass(r_in.copy(), 1, refpts=self.passrefs)
            trajs = r_out[:, 0, :, 0].T

        if not needs.isdisjoint({Need.ORBIT, Need.MATRIX, Need.LOCALOPTICS,
                                 Need.GLOBALOPTICS, Need.EMITTANCE}):
            o0, orbits = ring.find_orbit(refpts=self.orbitrefs)

        if Need.MATRIX in needs:
            if ring.is_6d:
                # noinspection PyUnboundLocalVariable
                _, mxdata = ring.find_m66(refpts=self.matrixrefs,
                                          orbit=o0, keep_lattice=True)
            else:
                # noinspection PyUnboundLocalVariable
                _, mxdata = ring.find_m44(refpts=self.matrixrefs,
                                          orbit=o0, keep_lattice=True)

        if not needs.isdisjoint({Need.LOCALOPTICS, Need.GLOBALOPTICS}):
            get_chrom = Need.CHROMATICITY in needs
            # noinspection PyUnboundLocalVariable
            _, rgdata, eldata = ring.get_optics(refpts=self.opticsrefs,
                                                orbit=o0, keep_lattice=True,
                                                get_chrom=get_chrom)

        if Need.EMITTANCE in needs:
            emdata = ring.envelope_parameters(orbit=o0, keep_lattice=True)

        for ob in self:
            eval(ob)

    # noinspection PyProtectedMember
    def exclude(self, obsname: str, excluded: Refpts):
        # Set the excluded mask on the selected observable
        for obs in self:
            if obs.name == obsname:
                obs._excluded = excluded
                obs._setup(self.ring)
        # Recompute all reference lists
        noref = self.ring.get_bool_index(None)
        self.orbitrefs = noref
        self.opticsrefs = noref
        self.passrefs = noref
        self.matrixrefs = noref
        for obs in self:
            self._update_reflists(obs)

    def clear(self):
        """Clear all the initial values"""
        for obs in self:
            obs.clear()

    def get_values(self, select: Optional[Container[str]] = None) -> list:
        """Return the values of selected observables

        Args:
            select:     :py:class:`~collections.abc.Container` of names for
              selecting observables. If :py:obj:`None` select all
        """
        selected = _selector(select)
        return [obs.value for obs in self if selected(obs)]

    values = property(get_values, doc="Values of the observables")

    def get_weights(self, select: Optional[Container[str]] = None) -> list:
        """Return the values of selected observables

        Args:
            select:     :py:class:`~collections.abc.Container` of names for
              selecting observables. If :py:obj:`None` select all
        """
        selected = _selector(select)
        return [np.broadcast_to(obs.w, np.asarray(obs.value).shape)
                for obs in self if selected(obs)]

    def get_weighted_values(self,
                            select: Optional[Container[str]] = None) -> list:
        """Return the values of selected observables

        Args:
            select:     :py:class:`~collections.abc.Container` of names for
              selecting observables. If :py:obj:`None` select all
        """
        selected = _selector(select)
        return [obs.weighted_value for obs in self if selected(obs)]

    def get_residuals(self, select: Optional[Container[str]] = None) -> list:
        """Return the residuals of selected observable

        Args:
            select:     :py:class:`~collections.abc.Container` of names for
              selecting observables. If :py:obj:`None` select all
        """
        selected = _selector(select)
        return [obs.residual for obs in self if selected(obs)]

    residuals = property(get_residuals, doc="Residuals of the observable")

    def get_flat_values(self,
                        select: Optional[Container[str]] = None) -> np.ndarray:
        """Return a 1-D array of selected Observable values

        Args:
            select:     :py:class:`~collections.abc.Container` of names for
              selecting observables. If :py:obj:`None` select all
        """
        selected = _selector(select)
        vals = (obs.value for obs in self if selected(obs))
        return np.concatenate([np.reshape(v, -1, order='F') for v in vals])

    flat_values = property(get_flat_values,
                           doc="1-D array of Observable values")

    def get_flat_weighted_values(self,
                                 select: Optional[Container[str]] = None
                                 ) -> np.ndarray:
        """Return a 1-D array of selected Observable weighted values

        Args:
            select:     :py:class:`~collections.abc.Container` of names for
              selecting observables. If :py:obj:`None` select all
        """
        selected = _selector(select)
        vals = (obs.weighted_value for obs in self if selected(obs))
        return np.concatenate([np.reshape(v, -1, order='F') for v in vals])

    flat_weighted_values = property(get_flat_weighted_values,
                                    doc="1-D array of Observable "
                                        "weigthed values")

    def get_flat_weights(self,
                         select: Optional[Container[str]] = None
                         ) -> np.ndarray:
        """Return a 1-D array of selected Observable weights

        Args:
            select:     :py:class:`~collections.abc.Container` of names for
              selecting observables. If :py:obj:`None` select all
        """
        vals = self.get_weights(select)
        return np.concatenate([np.reshape(v, -1, order='F') for v in vals])

    flat_weights = property(get_flat_weights,
                            doc="1-D array of Observable weights")

    def get_sum_residuals(self,
                          select: Optional[Container[str]] = None) -> float:
        """Return the sum of selected residual values

        Args:
            select:     :py:class:`~collections.abc.Container` of names for
              selecting observables. If :py:obj:`None` select all
        """
        selected = _selector(select)
        residuals = [obs.residual for obs in self if selected(obs)]
        return sum(np.sum(res) for res in residuals)

    sum_residuals = property(get_sum_residuals,
                             doc="Sum of all residual values")


# noinspection PyPep8Naming
def GlobalOpticsObservable(param: str, *args,
                           index: Optional[RefIndex] = Ellipsis,
                           name: Optional[str] = None,
                           use_integer: bool = False,
                           **kwargs):
    # noinspection PyUnresolvedReferences
    r"""Observe a global optics parameter

    Args:
        param:          Optics parameter name (see :py:func:`.get_optics`)
          or :ref:`user-defined evaluation function <globaloptics_eval>`
        index:          Index in the parameter array, If :py:obj:`Ellipsis`,
          the whole array is specified
        name:           Observable name. If :py:obj:`None`, an explicit
          name will be generated
        use_integer:    For *'tune'* parameter: derive the tune from the
          phase advance to avoid skipping integers. Slower than looking only
          at the fractional part
        *args:          Optional positional arguments for a user-defined
          evaluation function
        **kwargs:       Optional keyword arguments for a user-defined
          evaluation function

    Keyword Args:
        target:         Target value for a constraint. If :py:obj:`None`
          (default), the residual will always be zero.
        weight:         Weight factor: the residual is
          :pycode:`((value-target)/weight)**2`
        bounds:         Tuple of lower and upper bounds. The parameter
          is constrained in the interval
          [*target*\ -\ *low_bound* *target*\ +\ *up_bound*]

    The *target*, *weight* and *bounds* inputs must be broadcastable to the
    shape of *value*.

    .. _globaloptics_eval:
    .. rubric:: User-defined evaluation function

    It is called as:

    :pycode:`value = fun(ring, ringdata, *args, **kwargs)`

    *ringdata* if the output of :py:func:`.get_optics`,

    *value* is the value of the Observable.

    Examples:

        >>> obs = GlobalOpticsObservable('tune', use_integer=True)

        Observe the tune in both planes, including the integer part (slower)

        >>> obs = GlobalOpticsObservable('chromaticity', index=0)

        Observe the vertical chromaticity
    """
    if param == 'tune' and use_integer:
        def tune(ring, data):
            # noinspection PyProtectedMember
            mu = _ElementObservable._recordaccess('mu', index)(ring, data)
            return np.squeeze(mu, axis=0) / 2 / pi

        # noinspection PyProtectedMember
        name = _ElementObservable._set_name(name, tune, index)
        return LocalOpticsObservable(End, tune, *args,
                                     name=name,
                                     summary=True,
                                     use_integer=True, **kwargs)
    else:
        return _GlobalOpticsObservable(param, *args, index=index, name=name,
                                       **kwargs)
