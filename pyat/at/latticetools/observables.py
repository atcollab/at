r"""Definition of :py:class:`.Observable` objects used in matching and
response matrices.

Observables are ways to monitor various lattice parameters and use them in
matching and correction procedures.

Class hierarchy
---------------

:py:class:`Observable`\ (evalfun, name, target, weight, bounds, ...)
    :py:class:`RingObservable`\ (fun, ...)

    :py:func:`GlobalOpticsObservable`\ (attrname, plane, name, ...)

    :py:class:`EmittanceObservable`\ (attrname, plane, name, ...)

    :py:class:`ElementObservable`\ (fun, ...)
        :py:class:`OrbitObservable`\ (refpts, axis, name, ...)

        :py:class:`TrajectoryObservable`\ (refpts, axis, name, ...)

        :py:class:`MatrixObservable`\ (refpts, axis, name, ...)

        :py:class:`LocalOpticsObservable`\ (refpts, attrname, axis, name, ...)

        :py:class:`LatticeObservable`\ (refpts, attrname, index, name, ...)

        :py:class:`GeometryObservable`\ (refpts, attrname, name, ...)

:py:class:`.Observable`\ s are usually not evaluated directly, but through a
container which performs the required optics computation and feeds each
:py:class:`.Observable` with its specific data. After evaluation, each
:py:class:`.Observable` provides the following properties:

- :py:attr:`~Observable.value`
- :py:attr:`~Observable.weight`
- :py:attr:`~Observable.weighted_value`: :pycode:`value / weight`
- :py:attr:`~Observable.deviation`:  :pycode:`value - target`
- :py:attr:`~Observable.residual`:  :pycode:`((value - target)/weight)**2`
"""

from __future__ import annotations

__all__ = [
    "Observable",
    "RingObservable",
    "EmittanceObservable",
    "GlobalOpticsObservable",
    "ElementObservable",
    "OrbitObservable",
    "GeometryObservable",
    "LocalOpticsObservable",
    "LatticeObservable",
    "MatrixObservable",
    "TrajectoryObservable",
    "Need",
]

from collections.abc import Callable, Set
from enum import Enum
from itertools import repeat
from math import pi
from typing import Union

# For sys.version_info.minor < 9:
from typing import Tuple

import numpy as np

from ..lattice import AtError, AxisDef, axis_, plane_
from ..lattice import Lattice, Refpts, End

RefIndex = Union[int, Tuple[int, ...], slice]


# Observables must be pickleable. For this, the evaluation function must be a
# module-level function. No inner, nested function is allowed. So nested
# functions are replaced be module-level callable class instances:


class _ModFun:
    """General and pickleable evaluation function."""

    def __init__(self, fun, statfun):
        self.fun = fun
        self.statfun = statfun

    def __call__(self, *a):
        return self.statfun(self.fun(*a), axis=0)


class _ArrayAccess:
    """Access a selected item in an array."""

    def __init__(self, index):
        self.index = _all_rows(index)

    def __call__(self, data):
        index = self.index
        return data if index is None else data[self.index]


class _RecordAccess:
    """Access selected items in a record array."""

    def __init__(self, fieldname: str, index):
        self.index = index
        self.fieldname = fieldname

    def __call__(self, data):
        index = self.index
        data = getattr(data, self.fieldname)
        return data if index is None else data[self.index]


def _all_rows(index: RefIndex | None):
    """Prepends "all rows" (":") to an index tuple."""
    if index is None:
        return None
    if isinstance(index, tuple):
        return (slice(None),) + index
    else:
        return slice(None), index


class _Tune:
    """Get integer tune from the phase advance."""

    def __init__(self, idx: RefIndex):
        self.fun = _RecordAccess("mu", _all_rows(idx))

    def __call__(self, data):
        mu = self.fun(data)
        return np.squeeze(mu, axis=0) / 2 / pi


class _Ring:
    """Get an attribute of a lattice element."""

    def __init__(self, attrname, index, refpts):
        self.get_val = _RecordAccess(attrname, index)
        self.refpts = refpts

    def __call__(self, ring):
        vals = [self.get_val(el) for el in ring.select(self.refpts)]
        return np.array(vals)


class Need(Enum):
    """Defines the computation requirements for an :py:class:`Observable`."""

    #:  Provides the *ring* data to the evaluation function
    RING = 1
    #:  Specify :py:func:`.find_orbit` computation and provide its *orbit* output
    #:  to the evaluation function
    ORBIT = 2
    #:  Specify :py:func:`.find_m44` or :py:func:`.find_m44` computation and
    #:  provide its *m44* or *m66* output to the evaluation function
    MATRIX = 3
    #:  Specify :py:func:`.get_optics` computation and provide its *ringdata* output
    #:  to the evaluation function
    GLOBALOPTICS = 4
    #:  Specify :py:func:`.get_optics` computation and provide its *elemdata* output
    #:  to the evaluation function
    LOCALOPTICS = 5
    #:  Specify :py:func:`.lattice_pass` computation and provide its *r_out*
    #:  to the evaluation function
    TRAJECTORY = 6
    #:  Specify :py:func:`.envelope_parameters` computation and provide its
    #:  *params* output to the evaluation function
    EMITTANCE = 7
    #:  Associated with LOCALOPTICS, require local optics computation at all
    #:  points: slower but avoids jumps in phase advance
    ALL_POINTS = 8
    #:  Associated with LOCALOPTICS, require the *get_chrom* keyword
    CHROMATICITY = 9
    #:  Associatef with LOCALOPTICS, require the *get_w* keyword
    W_FUNCTIONS = 10
    #:  Specify :py:meth:`~.Lattice.get_geometry` computation and provide its ? output
    #:  to the evaluation function
    GEOMETRY = 11


class Observable:
    """Base class for Observables. Can be used for user-defined observables."""

    def __init__(
        self,
        fun: Callable,
        *args,
        name: str | None = None,
        target=None,
        weight=1.0,
        bounds=(0.0, 0.0),
        needs: Set[Need] | None = None,
        **kwargs,
    ):
        r"""Args:
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            fun:            :ref:`evaluation function <base_eval>`
            *args:          Arguments provided to the evaluation function
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

        Keyword Args:
            **kwargs:       Keyword arguments provided to the evaluation function

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        .. _base_eval:
        .. rubric:: User-defined evaluation function

        The general form is:

        :pycode:`value = fun(*data, *args, **kwargs)`

        - *data* depends on the *needs* argument, and by default is empty. If several
          needed values are specified, their order is: *ring*, *orbit*, *m44/m66*,
          *ringdata*, *elemdata*, *r_out*, *params*, *geomdata*,
        - *args* are the positional arguments provided to the observable constructor,
        - *kwargs* are the keyword arguments provided to the observable constructor,
        - *value* is the value of the observable.

        For user-defined evaluation functions using linear optics data or
        emittance data, it is recommended to use
        :py:class:`LocalOpticsObservable`, :py:obj:`GlobalOpticsObservable`
        or :py:class:`EmittanceObservable` which provide the corresponding
        *data* argument.
        """
        name = fun.__name__ if name is None else name
        self.fun: Callable = fun  #: Evaluation function
        self.needs: Set[Need] = needs or set()  #: Set of requirements
        self.name: str = name  #: Observable name
        self.target = target  #: Target value
        self.w = weight
        self.lbound, self.ubound = bounds
        self.initial = None
        self._value = None
        self._shape = None
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        """Return the string representation of the Observable."""
        return "\n".join((self._header(), self._all_lines()))

    @staticmethod
    def _header():
        """Header line."""
        fstring = "\n    {:<12} {:>16}  {:>16}  {:>16}  {:>16}  {:>16} "
        return fstring.format(
            "location", "Initial", "Actual", "Low bound", "High bound", "deviation"
        )

    @staticmethod
    def _line(loc, *items):
        def pitem(v):
            if v is None:
                return f" {'-':>14}   "
            elif isinstance(v, Exception):
                return f" {type(v).__name__:>16} "
            elif isinstance(v, np.ndarray) and v.ndim > 0:
                if v.size == 0:
                    return f" {'[]':16} "
                else:
                    return f" [{v.flat[0]:< 10.4} ...] "
            else:
                return f" {v: 16.6} "

        its = [pitem(v) for v in items]
        return f"    {loc:<12}" + "".join(its)

    def _all_lines(self):
        vnow = self._value
        if vnow is None or isinstance(vnow, Exception):
            deviation = None
        else:
            deviation = self.deviation
        if self.target is None:
            vmin = None
            vmax = None
        else:
            target = np.broadcast_to(self.target, np.asarray(vnow).shape)
            vmin = target + self.lbound
            vmax = target + self.ubound
        values = self._line("", self.initial, vnow, vmin, vmax, deviation)
        return "\n".join((self.name, values))

    def _setup(self, ring: Lattice):
        """Setup function called when the observable is added to a list."""
        pass

    def evaluate(self, *data, initial: bool = False):
        """Compute and store the value of the observable.

        The direct evaluation of a single :py:class:`Observable` is normally
        not used. This method is called by the :py:class:`.ObservableList`
        container which provides the *data* arguments.

        Args:
            *data:      Raw data, provided by :py:class:`.ObservableList` and
              sent to the evaluation function
            initial:    It :py:obj:`None`, store the result as the initial
              value
        """
        for d in data:
            if isinstance(d, Exception):
                self._value = d
                return d

        val = self.fun(*data, *self.args, **self.kwargs)
        if self._shape is None:
            self._shape = np.asarray(val).shape
        if initial:
            self.initial = val
        self._value = val
        return val

    @property
    def value(self):
        """Value of the observable."""
        val = self._value
        if isinstance(val, Exception):
            raise AtError(f"Evaluation of {self.name} failed: {val.args[0]}") from val
        return val

    @property
    def weight(self):
        """Observable weight."""
        return np.broadcast_to(self.w, np.asarray(self._value).shape)

    @property
    def weighted_value(self):
        """Weighted value of the Observable, computed as
        :pycode:`weighted_value = value/weight`.
        """
        return self.value / self.w

    @property
    def deviation(self):
        """Deviation from target value, computed as
        :pycode:`deviation = value-target`.
        """
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
        """:pycode:`weighted_deviation = (value-target)/weight`."""
        return self.deviation / self.w

    @property
    def residual(self):
        """residual, computed as :pycode:`residual = ((value-target)/weight)**2`."""
        return (self.deviation / self.w) ** 2

    @staticmethod
    def _set_name(name, param, index):
        """Compute a default observable names."""
        if name is None:
            if (
                index is Ellipsis
                or index is None
                or (isinstance(index, str) and index in {":", "..."})
            ):
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


class RingObservable(Observable):
    """Observe any user-defined property of a ring."""

    def __init__(
        self,
        fun: Callable,
        name: str | None = None,
        **kwargs,
    ):
        r"""Args:
            fun:            :ref:`user-defined evaluation function <ring_eval>`
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated.

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

        .. _ring_eval:
        .. rubric:: User-defined evaluation function

        It is called as:

        :pycode:`value = fun(ring)`

        - *ring* is the lattice description,
        - *value* is the value of the Observable.

        Examples:

            >>> def circumference(ring):
            ...     return ring.get_s_pos(len(ring))[0]
            >>> obs = RingObservable(circumference)

            Defines an Observable for the ring circumference.

            >>> def momentum_compaction(ring):
            ...     return ring.get_mcf()
            >>> obs = RingObservable(momentum_compaction)

            Defines an Observable for the momentum compaction factor.
        """
        needs = {Need.RING}
        name = fun.__name__ if name is None else name
        super().__init__(fun, name=name, needs=needs, **kwargs)


class ElementObservable(Observable):
    """Base class for Observables linked to a position in the lattice."""

    def __init__(
        self,
        fun: Callable,
        refpts: Refpts,
        name: str | None = None,
        summary: bool = False,
        statfun: Callable | None = None,
        **kwargs,
    ):
        r"""Args:
            fun:            :ref:`evaluation function <base_eval>`
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`.

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
            name = f"{statfun.__name__}({name})"
            fun = _ModFun(fun, statfun)
        super().__init__(fun, name=name, **kwargs)
        self.summary = summary
        self.refpts = refpts
        self._boolrefs = None
        self._excluded = None
        self._locations = [""]

    def _all_lines(self):
        if self.summary:
            return super()._all_lines()
        else:
            vnow = self._value
            if vnow is None or isinstance(vnow, Exception):
                vnow = repeat(vnow)
                deviation = repeat(None)
                vmin = repeat(None)
                vmax = repeat(None)
            else:
                deviation = self.deviation
                if self.target is None:
                    vmin = repeat(None)
                    vmax = repeat(None)
                else:
                    target = np.broadcast_to(self.target, np.asarray(vnow).shape)
                    vmin = target + self.lbound
                    vmax = target + self.ubound
            vini = self.initial
            if vini is None:
                vini = repeat(None)
            viter = zip(self._locations, vini, vnow, vmin, vmax, deviation)
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


class GeometryObservable(ElementObservable):
    """Observe the geometrical parameters of the reference trajectory.

    Process the *geomdata* output of :py:func:`.get_geometry`.
    """

    _field_list = {"x", "y", "angle"}

    def __init__(self, refpts: Refpts, param: str, name: str | None = None, **kwargs):
        # noinspection PyUnresolvedReferences
        r"""Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            param:          Geometry parameter name: one in {'x', 'y', 'angle'}
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated.

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

            >>> obs = GeometryObservable(at.Monitor, param="x")

            Observe x coordinate of monitors
        """
        if param not in self._field_list:
            raise ValueError(f"Expected {param!r} to be one of {self._field_list!r}")
        name = self._set_name(name, "geometry", param)
        fun = _RecordAccess(param, None)
        needs = {Need.GEOMETRY}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class OrbitObservable(ElementObservable):
    """Observe the closed orbit coordinates at selected locations.

    Process the *orbit* output of :py:func:`.find_orbit`.
    """

    def __init__(
        self, refpts: Refpts, axis: AxisDef = None, name: str | None = None, **kwargs
    ):
        # noinspection PyUnresolvedReferences
        r"""Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            axis:           Index in the orbit vector, If :py:obj:`None`,
              the whole vector is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated.

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
        name = self._set_name(name, "orbit", axis_(axis, "code"))
        fun = _ArrayAccess(axis_(axis, "index"))
        needs = {Need.ORBIT}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class MatrixObservable(ElementObservable):
    """Observe the closed orbit at selected locations.

    Processs the result of calling :py:func:`.find_m44` or :py:func:`.find_m44`
    depending upon :py:meth:`~.Lattice.is_6d`.
    """

    def __init__(
        self,
        refpts: Refpts,
        axis: AxisDef = Ellipsis,
        name: str | None = None,
        **kwargs,
    ):
        # noinspection PyUnresolvedReferences
        r"""Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            axis:           Index in the transfer matrix, If :py:obj:`Ellipsis`,
              the whole matrix is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated.

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

            >>> obs = MatrixObservable(at.Monitor, axis=("x", "px"))

            Observe the transfer matrix from origin to monitor locations and
            extract T[0,1]
        """
        name = self._set_name(name, "matrix", axis_(axis, "code"))
        fun = _ArrayAccess(axis_(axis, "index"))
        needs = {Need.MATRIX}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class _GlobalOpticsObservable(Observable):
    def __init__(
        self, param: str, plane: AxisDef = None, name: str | None = None, **kwargs
    ):
        # noinspection PyUnresolvedReferences
        r"""Args:
            param:          Optics parameter name (see :py:func:`.get_optics`)
              or user-defined evaluation function called as:
              :pycode:`value = fun(ringdata, ring=ring)` and returning the value of
              the Observable
            plane:          Index in the parameter array, If :py:obj:`None`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated.

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
        name = self._set_name(name, param, plane_(plane, "code"))
        if callable(param):
            fun = param
            needs.add(Need.CHROMATICITY)
        else:
            fun = _RecordAccess(param, plane_(plane, "index"))
            if param == "chromaticity":
                needs.add(Need.CHROMATICITY)
        super().__init__(fun, needs=needs, name=name, **kwargs)


class LocalOpticsObservable(ElementObservable):
    """Observe a local optics parameter at selected locations.

    Process the local output of :py:func:`.get_optics`.
    """

    def __init__(
        self,
        refpts: Refpts,
        param: str | Callable,
        plane: AxisDef = Ellipsis,
        name: str | None = None,
        all_points: bool = False,
        **kwargs,
    ):
        # noinspection PyUnresolvedReferences
        r"""Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            param:          Optics parameter name (see :py:func:`.get_optics`)
              or :ref:`user-defined evaluation function <localoptics_eval>`
            plane:          Index in the parameter array, If :py:obj:`Ellipsis`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            all_points:     Compute the local optics at all elements. This avoids
              discontinuities in phase advances. This is automatically set for the
              'mu' parameter, but may need to be specified for user-defined evaluation
              functions using the phase advance.

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

        :pycode:`value = fun(elemdata)`

        - *elemdata* is the output of :py:func:`.get_optics`, evaluated at the *refpts*
          of the observable,
        - *value* is the value of the Observable and must have one line per
          refpoint. Alternatively, it may be a single line, but then the
          *summary* keyword must be set to :py:obj:`True`.

        Examples:

            >>> obs = LocalOpticsObservable(at.Monitor, "beta")

            Observe the beta in both planes at all :py:class:`.Monitor`
            locations

            >>> obs = LocalOpticsObservable(
            ...     at.Quadrupole, "beta", plane="y", statfun=np.max
            ... )

            Observe the maximum vertical beta in Quadrupoles

            >>> def phase_advance(elemdata):
            ...     mu = elemdata.mu
            ...     return mu[-1] - mu[0]
            >>>
            >>> allobs.append(
            ...     LocalOpticsObservable(
            ...         [33, 101], phase_advance, all_points=True, summary=True
            ...     )
            ... )

            The user-defined evaluation function computes the phase-advance
            between the 1st and last given reference points, here the elements
            33 and 101 of the lattice
        """
        if param in {"M", "closed_orbit", "dispersion", "A", "R"}:
            ax_ = axis_
        else:
            ax_ = plane_

        needs = {Need.LOCALOPTICS}
        name = self._set_name(name, param, ax_(plane, "code"))
        if callable(param):
            fun = param
        else:
            fun = _RecordAccess(param, _all_rows(ax_(plane, "index")))
            if param == "mu" or all_points:
                needs.add(Need.ALL_POINTS)
            if param in {"W", "Wp", "dalpha", "dbeta", "dmu", "ddispersion", "dR"}:
                needs.add(Need.W_FUNCTIONS)

        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class LatticeObservable(ElementObservable):
    """Observe an attribute of selected lattice elements."""

    def __init__(
        self,
        refpts: Refpts,
        attrname: str,
        index: int | None = None,
        name: str | None = None,
        **kwargs,
    ):
        # noinspection PyUnresolvedReferences
        r"""Args:
            refpts:         Elements to be observed
              See ":ref:`Selecting elements in a lattice <refpts>`"
            attrname:       Attribute name
            index:          Index in the attribute array. If :py:obj:`None`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated.

        Keyword Args:
            statfun:        Post-processing function called on the value of the
              observable. Example: :pycode:`statfun=numpy.mean`

        Example:

            >>> obs = LatticeObservable(at.Sextupole, "KickAngle", index=0, statfun=np.sum)

            Observe the sum of horizontal kicks in Sextupoles
        """
        fun = _Ring(attrname, index, refpts)
        needs = {Need.RING}
        name = self._set_name(name, attrname, index)
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class TrajectoryObservable(ElementObservable):
    """Observe trajectory coordinates at selected locations.

    Process the *r_out* output if :py:meth:`.Lattice.track`
    """

    def __init__(
        self,
        refpts: Refpts,
        axis: AxisDef = Ellipsis,
        name: str | None = None,
        **kwargs,
    ):
        r"""Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            axis:          Index in the orbit array, If :py:obj:`Ellipsis`,
              the whole array is specified
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated.

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
        name = self._set_name(name, "trajectory", axis_(axis, "code"))
        fun = _ArrayAccess(axis_(axis, "index"))
        needs = {Need.TRAJECTORY}
        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)


class EmittanceObservable(Observable):
    """Observe emittance-related parameters.

    Process the output of :py:func:`.envelope_parameters`.
    """

    def __init__(
        self, param: str, plane: AxisDef = None, name: str | None = None, **kwargs
    ):
        r"""Args:
            param:          Parameter name (see :py:func:`.envelope_parameters`) or
              :ref:`user-defined evaluation function <emittance_eval>`
            plane:          One out of {0, 'x', 'h', 'H'} for horizontal plane,
             one out of {1, 'y', 'v', 'V'} for vertival plane or one out of
             {2, 'z', 'l', 'L'} for longitudinal plane
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated.

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

        .. _emittance_eval:
        .. rubric:: User-defined evaluation function

        It is called as:

        :pycode:`value = fun(paramdata)`

        *paramdata* if the :py:class:`.RingParameters` object returned by
        :py:func:`.envelope_parameters`.

        *value* is the value of the Observable.

        Example:

            >>> EmittanceObservable("emittances", plane="h")

            Observe the horizontal emittance
        """
        name = self._set_name(name, param, plane_(plane, "code"))
        if callable(param):
            fun = param
        else:
            fun = _RecordAccess(param, plane_(plane, "index"))
        needs = {Need.EMITTANCE}
        super().__init__(fun, needs=needs, name=name, **kwargs)


# noinspection PyPep8Naming
def GlobalOpticsObservable(
    param: str,
    plane: AxisDef = Ellipsis,
    name: str | None = None,
    use_integer: bool = False,
    **kwargs,
):
    # noinspection PyUnresolvedReferences
    r"""Observe a global optics parameter.

    Process the *ringdata* output of :py:func:`.get_optics`.

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

    - *ringdata* is the output of :py:func:`.get_optics`,
    - *value* is the value of the Observable.

    Examples:

        >>> obs = GlobalOpticsObservable("tune", use_integer=True)

        Observe the tune in both planes, including the integer part (slower)

        >>> obs = GlobalOpticsObservable("chromaticity", plane="v")

        Observe the vertical chromaticity
    """
    if param == "tune" and use_integer:
        # noinspection PyProtectedMember
        name = ElementObservable._set_name(name, param, plane_(plane, "code"))
        return LocalOpticsObservable(
            End,
            _Tune(plane_(plane, "index")),
            name=name,
            summary=True,
            all_points=True,
            **kwargs,
        )
    else:
        return _GlobalOpticsObservable(param, plane=plane, name=name, **kwargs)
