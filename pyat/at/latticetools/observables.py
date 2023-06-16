r"""
Definition of :py:class:`.Observable` objects used in matching and
response matrices

Observables are ways to monitor various lattice parameters and use them in
matching and correction procedures.

Class hierarchy
---------------

:py:class:`Observable`\ (evalfun, name, target, weight, bounds, ...)
    :py:func:`GlobalOpticsObservable`\ (attrname, plane, name, ...)

    :py:class:`EmittanceObservable`\ (attrname, plane, name, ...)

    :py:class:`_ElementObservable`\ (...)
        :py:class:`OrbitObservable`\ (refpts, axis, name, ...)

        :py:class:`BPMOrbitObservable`\ (refpts, axis, name, ...)

        :py:class:`TrajectoryObservable`\ (refpts, axis, name, ...)

        :py:class:`BPMTrajectoryObservable`\ (refpts, axis, name, ...)

        :py:class:`MatrixObservable`\ (refpts, axis, name, ...)

        :py:class:`LocalOpticsObservable`\ (refpts, attrname, axis, name, ...)

        :py:class:`LatticeObservable`\ (refpts, attrname, index, name, ...)

:py:class:`.Observable`\ s are usually not evaluated directly, but through a
container which performs the required optics computation and feeds each
:py:class:`.Observable` with its specific data. After evaluation, each
:py:class:`.Observable` provides the following properties:

- :py:attr:`~Observable.value`
- :py:attr:`~Observable.weight`
- :py:attr:`~Observable.weighted_value`: :pycode:`value / weight`
- :py:attr:`~Observable.deviation`:  :pycode:`value - target`
- :py:attr:`~Observable.residual`:  :pycode:`((value - target)/weight)**2`

:py:class:`ObservableList`\ (lattice, ...)
    This container based on :py:class:`list` is in charge of optics
    computations and provides each individual :py:class:`.Observable` with its
    specific data.

:py:class:`ObservableList` provides the :py:meth:`~ObservableList.evaluate`
method and the following properties:

- :py:attr:`~ObservableList.values`
- :py:attr:`~ObservableList.flat_values`
- :py:attr:`~ObservableList.weights`
- :py:attr:`~ObservableList.flat_weights`
- :py:attr:`~ObservableList.weighted_values`
- :py:attr:`~ObservableList.flat_weighted_values`
- :py:attr:`~ObservableList.deviations`
- :py:attr:`~ObservableList.flat_deviations`
- :py:attr:`~ObservableList.residuals`
- :py:attr:`~ObservableList.sum_residuals`

"""
from __future__ import annotations
from typing import Optional, Union
# For sys.version_info.minor < 9:
from typing import Tuple
from collections.abc import Callable, Iterable, Set
from functools import reduce
from math import pi
from enum import Enum
from itertools import repeat
from ..lattice import Lattice, Orbit, Refpts, All, End
from ..lattice import AxisDef, axis_, plane_, frequency_control
import numpy as np

RefIndex = Union[int, Tuple[int, ...], slice]

# Observables must be pickleable. For this, the evaluation function must be a
# module-level function. No inner, nested function is allowed. So nested
# functions are replaced be module-level callable class instances:


class _modfun(object):
    def __init__(self, fun, statfun):
        self.fun = fun
        self.statfun = statfun

    def __call__(self, ring, *a):
        return self.statfun(self.fun(ring, *a), axis=0)


class _arrayaccess(object):
    def __init__(self, index):
        self.index = _all_rows(index)

    def __call__(self, ring, data):
        index = self.index
        return data if index is None else data[self.index]


class _recordaccess(object):
    def __init__(self, fieldname, index):
        self.index = index
        self.fieldname = fieldname

    def __call__(self, ring, data):
        index = self.index
        data = getattr(data, self.fieldname)
        return data if index is None else data[self.index]


def _all_rows(index: RefIndex):
    if index is None:
        return None
    if isinstance(index, tuple):
        return (slice(None),) + index
    else:
        return slice(None), index


class _tune(object):
    def __init__(self, idx):
        self.fun = _recordaccess('mu', _all_rows(idx))

    def __call__(self, ring, data):
        mu = self.fun(ring, data)
        return np.squeeze(mu, axis=0) / 2 / pi


class _ring(object):
    def __init__(self, attrname, index, refpts):
        self.get_val = _recordaccess(attrname, index)
        self.refpts = refpts

    def __call__(self, ring):
        vals = [self.get_val(ring, el) for el in ring.select(self.refpts)]
        return np.array(vals)


def _flatten(vals, order='F'):
    return np.concatenate([np.reshape(v, -1, order=order) for v in vals])


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
    #:  to the evaluation function
    GLOBALOPTICS = 3
    #:  Specify :py:func:`.get_optics` computation and provide its *elemdata*
    #:  to the evaluation function
    LOCALOPTICS = 4
    #:  Specify :py:func:`.lattice_pass` computation and provide its *r_out*
    #:  to the evaluation function
    TRAJECTORY = 5
    #:  Specify :py:func:`.envelope_parameters` computation and provide its
    #:  *params* to the evaluation function
    EMITTANCE = 6
    #:  Associated with LOCALOPTICS, require local optics computation at all
    #:  points: slower but avoids jumps in phase advance
    ALL_POINTS = 7
    #:  Associated with LOCALOPTICS, require the *get_chrom* keyword
    CHROMATICITY = 8
    #:  Specify geometry computation and provide the full data at evaluation
    #:  points
    GEOMETRY = 9


class Observable(object):
    """Base class for Observables. Can be used for user-defined observables"""
    def __init__(self, fun: Callable,
                 name: Optional[str] = None,
                 target=None,
                 weight=1.0,
                 bounds=(0.0, 0.0),
                 needs: Optional[Set[Need]] = None,
                 fun_args: tuple = (),
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
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]
            needs:          Set of requirements. This selects the data provided
              to the evaluation function. *needs* items are members of the
              :py:class:`Need` enumeration

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        .. _base_eval:
        .. rubric:: User-defined evaluation function

        The general form is:

        :pycode:`value = fun(ring, *data)`

        *data* depends on the *needs* argument, and by default is empty.

        *value* is the value of the observable,

        For user-defined evaluation functions using linear optics data or
        emittance data, it is recommended to use
        :py:class:`LocalOpticsObservable`, :py:obj:`GlobalOpticsObservable`
        or :py:class:`EmittanceObservable` which provide the corresponding
        *data* argument.

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
        self.initial = None
        self._value = None
        self._shape = None
        self.args = fun_args
        self.kwargs = kwargs

    def __str__(self):
        return "\n".join((self._header(), self._all_lines()))

    @staticmethod
    def _header():
        """Header line"""
        fstring = "\n{:<12} {:>16}  {:>16}  {:>16}  {:>16}  {:>16} "
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
        return f"    {loc:<12}" + "".join(its)

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

        The direct evaluation of a single :py:class:`Observable` is normally
        not used. This method is called by the :py:class:`ObservableList`
        container which provides the *data* argument.

        Args:
            ring:       Lattice description
            *data:      Raw data, provided by :py:class:`ObservableList` and
              sent to the evaluation function
            initial:    It :py:obj:`None`, store the result as the initial
              value
        """
        val = self.fun(ring, *data, *self.fun_args, **self.kwargs)
        if self._shape is None:
            self._shape = np.asarray(val).shape
        if initial:
            self.initial = val
        self._value = val

    @property
    def value(self):
        """Value of the observable"""
        return self._value

    @property
    def weight(self):
        """Observable weight"""
        return np.broadcast_to(self.w, np.asarray(self.value).shape)

    @property
    def weighted_value(self):
        """Weighted value of the Observable, computed as
        :pycode:`deviation = value/weight`"""
        return self._value / self.w

    @property
    def deviation(self):
        """Deviation from target value, computed as
        :pycode:`deviation = value-target`
        """
        if self._value is None:
            deviation = None
        else:
            vnow = np.asarray(self.value)
            vsh = vnow.shape
            if self.target is None:
                deviation = np.broadcast_to(0.0, vsh)
            else:
                diff = np.atleast_1d(vnow - np.broadcast_to(self.target, vsh))
                lb = diff - self.lbound
                ub = diff - self.ubound
                lb[lb >= 0] = 0
                ub[ub <= 0] = 0
                deviation = (lb + ub).reshape(vsh)
        return deviation

    @property
    def weighted_deviation(self):
        return self.deviation / self.w

    @property
    def residual(self):
        """residual, computed as
        :pycode:`residual = ((value-target)/weight)**2`"""
        if self._value is None:
            res = None
        else:
            res = (self.deviation/self.w) ** 2
        return res

    @staticmethod
    def _set_name(name, param, index):
        """Compute a default observable names"""
        if name is None:
            if index is Ellipsis or index is None or \
                    (isinstance(index, str) and index in {":", "..."}):
                subscript = ""
            elif isinstance(index, tuple):
                ids = ", ".join(str(k) for k in index)
                subscript = f"[{ids}]"
            else:
                subscript = f"[{index}]"
            if callable(param):
                base = param.__name__
            else:
                base = param
            name = base + subscript
        return name


class _ElementObservable(Observable):
    """Base class for Observables linked to a position in the lattice"""
    def __init__(self, fun: Callable, refpts: Refpts,
                 name: Optional[str] = None,
                 summary: bool = False,
                 statfun: Optional[Callable] = None,
                 **kwargs):
        r"""
        Args:
            fun:            :ref:`evaluation function <base_eval>`
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`

        Keyword Args:
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.
        """
        name = fun.__name__ if name is None else name
        if statfun:
            summary = True
            name = '{}({})'.format(statfun.__name__, name)
            fun = _modfun(fun, statfun)
        super().__init__(fun, name=name, **kwargs)
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


class GeometryObservable(_ElementObservable):
    """Observe the geometrical parameters of the reference trajectory"""
    field_list = {'x', 'y', 'angle'}

    def __init__(self, refpts: Refpts, param: str,
                 name: Optional[str] = None,
                 **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            param:          Geometry parameter name: one in {'x', 'y', 'angle'}
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
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        Example:

            >>> obs = GeometryObservable(at.Monitor, param='x')

            Observe x coordinate of monitors
        """
        if param not in self.field_list:
            raise ValueError(
                f'Expected {param!r} to be one of {self.field_list!r}')
        name = self._set_name(name, 'geometry', param)
        fun = _recordaccess(param, None)
        needs = {Need.GEOMETRY}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class OrbitObservable(_ElementObservable):
    """Observes the transfer matrix at selected locations"""
    def __init__(self, refpts: Refpts, axis: AxisDef = None,
                 name: Optional[str] = None,
                 **kwargs):
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
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        Example:

            >>> obs = OrbitObservable(at.Monitor, axis=0)

            Observe the horizontal closed orbit at monitor locations
        """
        name = self._set_name(name, 'orbit', axis_(axis, 'code'))
        fun = _arrayaccess(axis_(axis, 'index'))
        needs = {Need.ORBIT}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class MatrixObservable(_ElementObservable):
    """Observes the closed orbit at selected locations"""
    def __init__(self, refpts: Refpts, axis: AxisDef = Ellipsis,
                 name: Optional[str] = None,
                 **kwargs):
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
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        Example:

            >>> obs = MatrixObservable(at.Monitor, axis=('x', 'px'))

            Observe the transfer matrix from origin to monitor locations and
            extract T[0,1]
        """
        name = self._set_name(name, 'matrix', axis_(axis, 'code'))
        fun = _arrayaccess(axis_(axis, 'index'))
        needs = {Need.MATRIX}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class _GlobalOpticsObservable(Observable):
    def __init__(self, param: str, plane: AxisDef = None,
                 name: Optional[str] = None,
                 **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            param:          Optics parameter name (see :py:func:`.get_optics`)
              or user-defined evaluation function called as:
              :pycode:`value = fun(ring, ringdata)` and returning the value of
              the Observable
            plane:          Index in the parameter array, If :py:obj:`None`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated

        Keyword Args:
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.
        """
        needs = {Need.GLOBALOPTICS}
        name = self._set_name(name, param, plane_(plane, 'code'))
        if callable(param):
            fun = param
            needs.add(Need.CHROMATICITY)
        else:
            fun = _recordaccess(param, plane_(plane, 'index'))
            if param == 'chromaticity':
                needs.add(Need.CHROMATICITY)
        super().__init__(fun, needs=needs, name=name, **kwargs)


class LocalOpticsObservable(_ElementObservable):
    """Observe a local optics parameter at selected locations"""
    def __init__(self, refpts: Refpts, param: Union[str, Callable],
                 plane: AxisDef = Ellipsis,
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
            plane:          Index in the parameter array, If :py:obj:`Ellipsis`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            use_integer:    For  the *'mu'* parameter, compute the
              phase advance at all points to avoid discontinuities (slower)

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
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        .. _localoptics_eval:
        .. rubric:: User-defined evaluation function

        It is called as:

        :pycode:`value = fun(ring, elemdata)`

        *elemdata* if the output of :py:func:`.get_optics`.

        *value* is the value of the Observable and must have one line per
        refpoint. Alternatively, it may be a single line, but then the
        *summary* keyword must be set to :py:obj:`True`

        Examples:

            >>> obs = LocalOpticsObservable(at.Monitor, 'beta')

            Observe the beta in both planes at all :py:class:`.Monitor`
            locations

            >>> obs = LocalOpticsObservable(at.Quadrupole, 'beta', plane='y',
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
        if param in {'M', 'closed_orbit', 'dispersion', 'A', 'R'}:
            ax_ = axis_
        else:
            ax_ = plane_
        needs = {Need.LOCALOPTICS}
        name = self._set_name(name, param, ax_(plane, 'code'))
        if callable(param):
            fun = param
            needs.add(Need.CHROMATICITY)
        else:
            fun = _recordaccess(param, _all_rows(ax_(plane, 'index')))
        if use_integer:
            needs.add(Need.ALL_POINTS)

        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class LatticeObservable(_ElementObservable):
    """Observe an attribute of selected lattice elements"""
    def __init__(self, refpts: Refpts, attrname: str, index: AxisDef = Ellipsis,
                 name: Optional[str] = None,
                 **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            refpts:         Elements to be observed
              See ":ref:`Selecting elements in a lattice <refpts>`"
            attrname:       Attribute name
            index:          Index in the attribute array. If :py:obj:`Ellipsis`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated

        Keyword Args:
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`

        Example:

            >>> obs = LatticeObservable(at.Sextupole, 'KickAngle', plane=0,
            ...                      statfun=np.sum)

            Observe the sum of horizontal kicks in Sextupoles
        """
        fun = _ring(attrname, index, refpts)
        name = self._set_name(name, attrname, index)
        super().__init__(fun, refpts, name=name, **kwargs)


class TrajectoryObservable(_ElementObservable):
    """Observe trajectory coordinates at selected locations"""
    def __init__(self, refpts: Refpts, axis: AxisDef = Ellipsis,
                 name: Optional[str] = None,
                 **kwargs):
        r"""
        Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            axis:          Index in the orbit array, If :py:obj:`Ellipsis`,
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
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.
        """
        name = self._set_name(name, 'trajectory', axis_(axis, 'code'))
        fun = _arrayaccess(axis_(axis, 'index'))
        needs = {Need.TRAJECTORY}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class EmittanceObservable(Observable):
    """Observe emittance-related parameters"""
    def __init__(self, param: str, plane: AxisDef = None,
                 name: Optional[str] = None,
                 **kwargs):
        r"""
        Args:
            param:          Parameter name (see
              :py:func:`.envelope_parameters`)
            plane:          One out of {0, 'x', 'h', 'H'} for horizontal plane,
             one out of {1, 'y', 'v', 'V'} for vertival plane or one out of
             {2, 'z', 'l', 'L'} for longitudinal plane
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
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        Example:

            >>> EmittanceObservable('emittances', plane='h')

            Observe the horizontal emittance
        """
        name = self._set_name(name, param, plane_(plane, 'code'))
        fun = _recordaccess(param, plane_(plane, 'index'))
        needs = {Need.EMITTANCE}
        super().__init__(fun, needs=needs, name=name, **kwargs)


class ObservableList(list):
    """Handles a list of Observables to be evaluated together

    :py:class:`ObservableList` supports all :py:class:`list` methods, like
    appending, insertion or concatenation with the "+" operator.
    """
    def __init__(self, ring: Lattice, obsiter: Iterable[Observable] = (),
                 **kwargs):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            ring:       Lattice description. This lattice is used for sorting
              the locations of the different observables but is not used for
              evaluation: the lattice used for evaluation is given as argument
              to the :py:meth:`~ObservableList.evaluate` method. The design
              lattice may be used here, for instance.
            obsiter:    Iterable of :py:class:`Observable`\ s

        Example:

            >>> obslist = ObservableList(ring)

            Create an empty Observable list

            >>> obslist.append(OrbitObservable(at.Monitor, plane='x'))
            >>> obslist.append(GlobalOpticsObservable('tune'))
            >>> obslist.append(EmittanceObservable('emittances', plane='h'))

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

            >>> obslist.get_flat_values({'tune', 'emittances[h]'})
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
        self.kwargs = kwargs
        super().__init__(self._scan(obs) for obs in obsiter)

    # noinspection PyProtectedMember
    def _update_reflists(self, obs) -> None:
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

    def _scan(self, obs) -> Observable:
        # When unpickling the ObservableList, the list is built before
        # self.__dict__ is restored, so the "ring" attribute is missing.
        # The "needs" and reflists will be restored later
        if hasattr(self, 'ring'):
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

    def __add__(self, other) -> ObservableList:
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

    def evaluate(self, ring: Lattice, *,
                 r_in: Orbit = None,
                 initial: bool = False, **kwargs):
        r"""Compute all the :py:class:`Observable` values

        Args:
            ring:       Lattice description used for evaluation
            r_in:       Optional coordinate input for trajectory observables
            initial:    If :py:obj:`True`, store the values as *initial values*

        Keyword Args:
            dp (float):     Momentum deviation. Defaults to :py:obj:`None`
            dct (float):    Path lengthening. Defaults to :py:obj:`None`
            df (float):     Deviation from the nominal RF frequency.
              Defaults to :py:obj:`None`
        """
        def obseval(obs):
            """Evaluate a single observable"""
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
            if Need.GEOMETRY in obsneeds:
                data.append(geodata[obsrefs])
            obs.evaluate(ring, *data, initial=initial)

        @frequency_control
        def ringeval(ring, dp: Optional[float] = None,
                     dct: Optional[float] = None,
                     df: Optional[float] = None,
                     r_in: Orbit = None):
            """Optics computations"""

            trajs = orbits = rgdata = eldata = emdata = mxdata = geodata = None
            needs = self.needs

            if Need.TRAJECTORY in needs:
                if r_in is None:
                    r_in = np.zeros(6)
                r_out = ring.lattice_pass(r_in.copy(), 1, refpts=self.passrefs)
                trajs = r_out[:, 0, :, 0].T

            if not needs.isdisjoint({Need.ORBIT, Need.MATRIX, Need.LOCALOPTICS,
                                     Need.GLOBALOPTICS, Need.EMITTANCE}):
                orbit0 = self.kwargs.pop('orbit', None)
                o0, orbits = ring.find_orbit(refpts=self.orbitrefs,
                                             dp=dp, dct=dct, df=df,
                                             orbit=orbit0)

            if Need.MATRIX in needs:
                if ring.is_6d:
                    # noinspection PyUnboundLocalVariable
                    _, mxdata = ring.find_m66(refpts=self.matrixrefs,
                                              dp=dp, dct=dct, df=df,
                                              orbit=o0, keep_lattice=True)
                else:
                    # noinspection PyUnboundLocalVariable
                    _, mxdata = ring.find_m44(refpts=self.matrixrefs,
                                              dp=dp, dct=dct, df=df,
                                              orbit=o0, keep_lattice=True)

            if not needs.isdisjoint({Need.LOCALOPTICS, Need.GLOBALOPTICS}):
                get_chrom = Need.CHROMATICITY in needs
                twiss_in = self.kwargs.pop('twiss_in', None)
                # noinspection PyUnboundLocalVariable
                _, rgdata, eldata = ring.get_optics(refpts=self.opticsrefs,
                                                    dp=dp, dct=dct, df=df,
                                                    orbit=o0, keep_lattice=True,
                                                    get_chrom=get_chrom,
                                                    twiss_in=None)

            if Need.EMITTANCE in needs:
                emdata = ring.envelope_parameters(orbit=o0, keep_lattice=True)

            if Need.GEOMETRY in needs:
                geodata, _ = ring.get_geometry()

            return trajs, orbits, rgdata, eldata, emdata, mxdata, geodata

        trajs, orbits, rgdata, eldata, emdata, mxdata, geodata = \
            ringeval(ring, r_in=r_in, **kwargs)

        for ob in self:
            obseval(ob)

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

    @property
    def shapes(self) -> list:
        """Shapes of all values"""
        # noinspection PyProtectedMember
        return [obs._shape for obs in self]

    @property
    def flat_shape(self):
        """Shape of the flattened values"""
        # noinspection PyProtectedMember
        vals = (reduce(lambda x, y: x*y, obs._shape, 1) for obs in self)
        return sum(vals),

    @property
    def values(self) -> list:
        """Return the values of all observables
        """
        return [obs.value for obs in self]

    def get_flat_values(self, order: str = 'F') -> np.ndarray:
        """Return a 1-D array of all Observable values

        Args:
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten((obs.value for obs in self), order=order)

    @property
    def weighted_values(self) -> list:
        """Weighted values of all observables"""
        return [obs.weighted_value for obs in self]

    def get_flat_weighted_values(self, order: str = 'F') -> np.ndarray:
        """Return a 1-D array of all Observable weighted values

        Args:
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten((obs.weighted_value for obs in self), order=order)

    @property
    def deviations(self) -> list:
        """Deviations from target values"""
        return [obs.deviation for obs in self]

    def get_flat_deviations(self, order: str = 'F') -> np.ndarray:
        """Return a 1-D array of deviations from target values

        Args:
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten((obs.deviation for obs in self), order=order)

    @property
    def weighted_deviations(self) -> list:
        """Weighted deviations from target values"""
        return [obs.weighted_deviation for obs in self]

    def get_flat_weighted_deviations(self, order: str = 'F') -> np.ndarray:
        """Return a 1-D array of weighted deviations from target values

        Args:
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten((obs.weighted_deviation for obs in self), order=order)

    @property
    def weights(self) -> list:
        """Weights of all observables"""
        return [obs.weight for obs in self]

    def get_flat_weights(self, order: str = 'F') -> np.ndarray:
        """Return a 1-D array of all Observable weights

        Args:
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten((obs.weight for obs in self), order=order)

    @property
    def residuals(self) -> list:
        """Residuals of all observable"""
        return [obs.residual for obs in self]

    def get_sum_residuals(self) -> float:
        """Return the sum of all residual values
        """
        residuals = (obs.residual for obs in self)
        return sum(np.sum(res) for res in residuals)

    flat_values = property(get_flat_values,
                           doc="1-D array of Observable values")
    flat_weighted_values = property(get_flat_weighted_values,
                                    doc="1-D array of Observable "
                                        "weigthed values")
    flat_deviations = property(get_flat_deviations,
                               doc="1-D array of deviations from target value")
    flat_weighted_deviations = property(get_flat_weighted_deviations)
    flat_weights = property(get_flat_weights,
                            doc="1-D array of Observable weights")
    sum_residuals = property(get_sum_residuals,
                             doc="Sum of all residual values")


# noinspection PyPep8Naming
def GlobalOpticsObservable(param: str, plane: AxisDef = Ellipsis,
                           name: Optional[str] = None,
                           use_integer: bool = False,
                           **kwargs):
    # noinspection PyUnresolvedReferences
    r"""Observe a global optics parameter

    Args:
        param:          Optics parameter name (see :py:func:`.get_optics`)
          or :ref:`user-defined evaluation function <globaloptics_eval>`
        plane:          Index in the parameter array, If :py:obj:`Ellipsis`,
          the whole array is specified
        name:           Observable name. If :py:obj:`None`, an explicit
          name will be generated
        use_integer:    For *'tune'* parameter: derive the tune from the
          phase advance to avoid skipping integers. Slower than looking only
          at the fractional part

    Keyword Args:
        target:         Target value for a constraint. If :py:obj:`None`
          (default), the residual will always be zero.
        weight:         Weight factor: the residual is
          :pycode:`((value-target)/weight)**2`
        bounds:         Tuple of lower and upper bounds. The parameter
          is constrained in the interval
          [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

    The *target*, *weight* and *bounds* inputs must be broadcastable to the
    shape of *value*.

    .. _globaloptics_eval:
    .. rubric:: User-defined evaluation function

    It is called as:

    :pycode:`value = fun(ring, ringdata)`

    *ringdata* if the output of :py:func:`.get_optics`,

    *value* is the value of the Observable.

    Examples:

        >>> obs = GlobalOpticsObservable('tune', use_integer=True)

        Observe the tune in both planes, including the integer part (slower)

        >>> obs = GlobalOpticsObservable('chromaticity', plane='v')

        Observe the vertical chromaticity
    """
    if param == 'tune' and use_integer:
        # noinspection PyProtectedMember
        name = _ElementObservable._set_name(name, param, plane_(plane, 'code'))
        return LocalOpticsObservable(End, _tune(plane_(plane, 'index')),
                                     name=name, summary=True, use_integer=True,
                                     **kwargs)
    else:
        return _GlobalOpticsObservable(param, plane=plane, name=name, **kwargs)
