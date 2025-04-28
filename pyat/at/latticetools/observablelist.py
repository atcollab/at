r"""Grouping of :py:class:`.Observable` objects for fast evaluation.

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

__all__ = [
    "ObservableList",
]

from collections.abc import Iterable, Iterator
from functools import reduce
from typing import Callable

import numpy as np
import numpy.typing as npt

# noinspection PyProtectedMember
from .observables import Observable, ElementObservable, Need
from .rdt_observable import RDTObservable
from ..lattice import AtError, frequency_control
from ..lattice import Lattice, Orbit, Refpts, All
from ..physics import linopt6
from ..tracking import internal_lpass


def _flatten(vals, order="F") -> npt.NDArray[float]:
    return np.concatenate([np.reshape(v, -1, order=order) for v in vals])


class _ObsResIter(Iterator):
    """Iterator object for the _ObsResult class"""

    def __init__(self, obsiter):
        self.base = obsiter

    def __next__(self):
        # Raises the stored error when reaching a missing value
        return Observable.check_value(next(self.base))


class _ObsResults(tuple):
    """Tuple-like object for the output of  ObservableList.evaluate

    _ObsResult implements a special treatment when the evaluation ends with an error.
    The error is stored instead of the value. This object raises the error a posteriori,
    when accessing the missing value.
    """

    def __getitem__(self, item):
        # Raises the stored error when accessing a missing value
        if isinstance(item, slice):
            return _ObsResults(super().__getitem__(item))
        else:
            return Observable.check_value(super().__getitem__(item))

    def __iter__(self):
        return _ObsResIter(super().__iter__())

    def __repr__(self):
        return repr(tuple(super().__iter__()))


class ObservableList(list):
    """Handles a list of Observables to be evaluated together.

    :py:class:`ObservableList` supports all :py:class:`list` methods, like
    appending, insertion or concatenation with the "+" operator.
    """

    needs_ring = {
        Need.RING,
        Need.ORBIT,
        Need.MATRIX,
        Need.GLOBALOPTICS,
        Need.LOCALOPTICS,
        Need.TRAJECTORY,
        Need.EMITTANCE,
        Need.GEOMETRY,
    }
    needs_orbit = {
        Need.ORBIT,
        Need.MATRIX,
        Need.GLOBALOPTICS,
        Need.LOCALOPTICS,
        Need.EMITTANCE,
    }
    needs_optics = {Need.GLOBALOPTICS, Need.LOCALOPTICS}

    def __init__(
        self,
        obsiter: Iterable[Observable] = (),
        *,
        method: Callable = linopt6,
        orbit: Orbit = None,
        twiss_in=None,
        r_in: Orbit = None,
        **kwargs,
    ):
        # noinspection PyUnresolvedReferences
        r"""
        Args:
            obsiter:    Iterable of :py:class:`.Observable`\ s

        Keyword Args:
            orbit (Orbit):      Initial orbit. Avoids looking for the closed
              orbit if it is already known. Used for
              :py:class:`MatrixObservable` and :py:class:`LocalOpticsObservable`
            twiss_in:           Initial conditions for transfer line optics.
              See :py:func:`.get_optics`. Used for
              :py:class:`LocalOpticsObservable`
            method (Callable):  Method for linear optics. Used for
              :py:class:`LocalOpticsObservable`.
              Default: :py:obj:`~.linear.linopt6`
            r_in (Orbit):       Initial trajectory, used for
              :py:class:`TrajectoryObservable`, Default: zeros(6)


        Example:
            >>> obslist = ObservableList()

            Create an empty Observable list

            >>> obslist.append(OrbitObservable(at.Monitor, plane="x"))
            >>> obslist.append(GlobalOpticsObservable("tune"))
            >>> obslist.append(EmittanceObservable("emittances", plane="h"))

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

            >>> obslist.get_flat_values("tune", "emittances[h]")
            array([3.815630e-01, 8.543764e-01, 1.090607e-04, 1.320391e-10])

            Get a flattened array of tunes and horizontal emittance
        """
        self.orbitrefs = None
        self.opticsrefs = None
        self.passrefs = None
        self.matrixrefs = None
        self.rdtrefs = None
        self.needs = None
        self.rdt_type = set()
        self.method = method
        self.orbit = orbit
        self.twiss_in = twiss_in
        self.r_in = r_in
        self.kwargs = kwargs
        super().__init__(obsiter)

    # noinspection PyProtectedMember
    def _setup(self, ring: Lattice):
        # Compute the union of all needs
        needs = set()
        rdt_type = set()
        for obs in self:
            needs |= obs.needs
            if isinstance(obs, ElementObservable):
                needs.add(Need.RING)
            if isinstance(obs, RDTObservable):
                rdt_type.add(obs._rdt_type)

        if (needs & self.needs_ring) and ring is None:
            raise ValueError("At least one Observable needs a ring argument")

        self.needs = needs
        self.rdt_type = rdt_type
        if ring is None:
            # Initialise each observable
            for obs in self:
                obs._setup(ring)
        else:
            # Initialise each observable and make a summary all refpoints
            noref = ring.get_bool_index(None)
            orbitrefs = opticsrefs = passrefs = matrixrefs = rdtrefs = noref
            for obs in self:
                obs._setup(ring)
                obsneeds = obs.needs
                if isinstance(obs, ElementObservable):
                    if Need.ORBIT in obsneeds:
                        orbitrefs |= obs._boolrefs
                    if Need.MATRIX in obsneeds:
                        matrixrefs |= obs._boolrefs
                    if Need.LOCALOPTICS in obsneeds:
                        if Need.ALL_POINTS in obsneeds:
                            opticsrefs = ring.get_bool_index(All)
                        else:
                            opticsrefs |= obs._boolrefs
                    if Need.TRAJECTORY in obsneeds:
                        passrefs |= obs._boolrefs
                    if Need.RDT in obsneeds:
                        rdtrefs |= obs._boolrefs
            self.orbitrefs = orbitrefs
            self.opticsrefs = opticsrefs
            self.rdtrefs = rdtrefs
            self.passrefs = passrefs
            self.matrixrefs = matrixrefs

    def __iadd__(self, other: ObservableList):
        if not isinstance(other, ObservableList):
            raise TypeError(f"Cannot add a {type(other)} to an ObservableList")
        self.extend(other)
        return self

    def __add__(self, other) -> ObservableList:
        nobs = ObservableList(self)
        nobs += other
        return nobs

    def append(self, obs: Observable):
        """Append observable to the end of the list."""
        if not isinstance(obs, Observable):
            raise TypeError(f"Cannot append a {type(obs)} to an ObservableList")
        self.needs = None
        super().append(obs)

    def extend(self, obsiter: Iterable[Observable]):
        """Extend list by appending Observables from the iterable."""
        self.needs = None
        super().extend(obsiter)

    def insert(self, index: int, obs: Observable):
        """Insert Observable before index."""
        if not isinstance(obs, Observable):
            raise TypeError(f"Cannot insert a {type(obs)} in an ObservableList")
        self.needs = None
        super().insert(index, obs)

    # noinspection PyProtectedMember
    def __str__(self):
        values = "\n".join(obs._all_lines() for obs in self)
        return "\n".join((Observable._header(), values))

    def evaluate(
        self,
        ring: Lattice | None = None,
        *,
        dp: float | None = None,
        dct: float | None = None,
        df: float | None = None,
        initial: bool = False,
        **kwargs,
    ):
        r"""Compute all the :py:class:`Observable` values.

        Args:
            ring:           Lattice used for evaluation
            dp (float):     Momentum deviation. Defaults to :py:obj:`None`
            dct (float):    Path lengthening. Defaults to :py:obj:`None`
            df (float):     Deviation from the nominal RF frequency.
              Defaults to :py:obj:`None`
            initial:    If :py:obj:`True`, store the values as *initial values*

        Keyword Args:
            orbit (Orbit):  Initial orbit. Avoids looking for the closed
              orbit if it is already known. Used for
              :py:class:`.MatrixObservable` and :py:class:`.LocalOpticsObservable`
            twiss_in:       Initial conditions for transfer line optics.
              See :py:func:`.get_optics`. Used for
              :py:class:`.LocalOpticsObservable`
            method (Callable):  Method for linear optics. Used for
              :py:class:`.LocalOpticsObservable`.
              Default: :py:obj:`~.linear.linopt6`
            r_in (Orbit):   Initial trajectory, used for
              :py:class:`.TrajectoryObservable`, Default: zeros(6)
        """

        def obseval(ring, obs):
            """Evaluate a single observable."""

            def check_error(data, refpts):
                return data if isinstance(data, Exception) else data[refpts]

            obsneeds = obs.needs
            obsrefs = getattr(obs, "_boolrefs", None)
            data = []
            if Need.RING in obsneeds:
                data.append(ring)
            if Need.ORBIT in obsneeds:
                data.append(check_error(orbits, obsrefs[self.orbitrefs]))
            if Need.MATRIX in obsneeds:
                data.append(check_error(mxdata, obsrefs[self.matrixrefs]))
            if Need.GLOBALOPTICS in obsneeds:
                data.append(rgdata)
            if Need.LOCALOPTICS in obsneeds:
                data.append(check_error(eldata, obsrefs[self.opticsrefs]))
            if Need.TRAJECTORY in obsneeds:
                data.append(trajs[obsrefs[self.passrefs]])
            if Need.EMITTANCE in obsneeds:
                data.append(emdata)
            if Need.GEOMETRY in obsneeds:
                data.append(geodata[obsrefs])
            if Need.RDT in obsneeds:
                data.append(check_error(rdtdata, obsrefs[self.rdtrefs]))
            return obs.evaluate(*data, initial=initial)

        @frequency_control
        def ringeval(
            ring,
            dp: float | None = None,
            dct: float | None = None,
            df: float | None = None,
        ):
            """Optics computations."""
            keep_lattice = False
            trajs = orbits = rgdata = eldata = emdata = mxdata = geodata = rdtdata = (
                None
            )
            twiss_in = kwargs.get("twiss_in", self.twiss_in)
            o0 = kwargs.get("orbit", self.orbit)
            o0 = getattr(twiss_in, "closed_orbit", None) if o0 is None else o0
            needs = self.needs
            needs_o0 = (needs & self.needs_orbit) and (o0 is None)

            if Need.TRAJECTORY in needs:
                # Trajectory computation
                r_in = kwargs.get("r_in", self.r_in)
                if r_in is None:
                    r_in = np.zeros(6)
                r_out = internal_lpass(ring, r_in.copy(), 1, refpts=self.passrefs)
                trajs = r_out[:, 0, :, 0].T
                keep_lattice = True

            # if needs & self.needs_orbit:
            if Need.ORBIT in needs or needs_o0:
                # Closed orbit computation
                try:
                    o0, orbits = ring.find_orbit(
                        refpts=self.orbitrefs,
                        dp=dp,
                        dct=dct,
                        df=df,
                        orbit=o0,
                        keep_lattice=keep_lattice,
                    )
                except AtError as err:
                    orbits = mxdata = rgdata = eldata = emdata = err
                else:
                    keep_lattice = True

            if Need.MATRIX in needs and o0 is not None:
                # Transfer matrix computation
                find_m = ring.find_m66 if ring.is_6d else ring.find_m44
                # noinspection PyUnboundLocalVariable
                _, mxdata = find_m(
                    refpts=self.matrixrefs,
                    dp=dp,
                    dct=dct,
                    df=df,
                    orbit=o0,
                    keep_lattice=keep_lattice,
                )
                keep_lattice = True

            if (needs & self.needs_optics) and o0 is not None:
                # Linear optics computation
                try:
                    _, rgdata, eldata = ring.get_optics(
                        refpts=self.opticsrefs,
                        dp=dp,
                        dct=dct,
                        df=df,
                        orbit=o0,
                        keep_lattice=keep_lattice,
                        get_chrom=Need.CHROMATICITY in needs,
                        get_w=Need.W_FUNCTIONS in needs,
                        twiss_in=twiss_in,
                        method=kwargs.get("method", self.method),
                    )
                except AtError as err:
                    rgdata = eldata = err
                else:
                    keep_lattice = True

            if Need.EMITTANCE in needs and o0 is not None:
                # Emittance computation
                try:
                    emdata = ring.envelope_parameters(
                        orbit=o0, keep_lattice=keep_lattice
                    )
                except Exception as err:
                    emdata = err

            if Need.GEOMETRY in needs:
                # Geometry computation
                geodata, _ = ring.get_geometry()

            if Need.RDT in needs:
                # RDT computation
                use_mp = kwargs.get("use_mp", False)
                pool_size = kwargs.get("pool_size", None)
                try:
                    _, _, rdtdata = ring.get_rdts(
                        refpts=self.rdtrefs,
                        rdt_type=self.rdt_type,
                        second_order=Need.RDT_2ND_ORDER in needs,
                        use_mp=use_mp,
                        pool_size=pool_size,
                    )
                except Exception as err:
                    rdtdata = err

            return trajs, orbits, rgdata, eldata, emdata, mxdata, geodata, rdtdata

        if self.needs is None or initial:
            self._setup(ring)

        trajs, orbits, rgdata, eldata, emdata, mxdata, geodata, rdtdata = ringeval(
            ring, dp=dp, dct=dct, df=df
        )
        return _ObsResults(obseval(ring, ob) for ob in self)

    def check(self) -> bool:
        """Check if all observables are evaluated.

        Returns:
            ok: :py:obj:`True` if evaluation is done, :py:obj:`False` otherwise

        Raises:
            AtError:    any value is doubtful: evaluation failed, empty value…
        """
        return all(obs.check() for obs in self)

    # noinspection PyProtectedMember
    def exclude(self, obsname: str, excluded: Refpts):
        """Set the excluded mask on the selected observable."""
        for obs in self:
            if obs.name == obsname:
                obs._excluded = excluded
        self.needs = None

    def _lookup(self, *ids: int | str) -> list[Observable]:
        """Observable lookup function"""

        def select(id):
            if isinstance(id, str):
                for obs in self:
                    if obs.name == id:
                        return obs
                else:
                    raise KeyError(id)
            else:
                return self[id]

        if ids:
            return [select(id) for id in ids]
        else:
            return self

    def _collect(self, attrname: str, *obsid: str | int, err: float | None = None):
        def val(obs):
            try:
                res = getattr(obs, attrname)
            except AtError:
                if err is None:
                    raise
                else:
                    # noinspection PyProtectedMember
                    shp = obs._shape
                    res = err if shp is None else np.broadcast_to(err, shp)
            return res

        obslist = self._lookup(*obsid)
        return tuple(val(obs) for obs in obslist)

    def get_shapes(self, *obsid: str | int) -> tuple:
        """Return the shapes of all values.

        Args:
            *obsid: name or index of selected observables (Default all)
        """
        return self._collect("_shape", *obsid)

    def get_flat_shape(self, *obsid: str | int):
        """Return the shape of the flattened values.

        Args:
            *obsid: name or index of selected observables (Default all)
        """
        vals = (
            reduce(lambda x, y: x * y, shp, 1)
            for shp in self._collect("_shape", *obsid)
        )
        return (sum(vals),)

    def get_values(self, *obsid: str | int, err: float | None = None) -> tuple:
        """Return the values of observables.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.

        Raises:
            Exception: Any exception raised during evaluation, unless *err* has been
              set.
        """
        return self._collect("value", *obsid, err=err)

    def get_flat_values(
        self, *obsid: str | int, err: float | None = None, order: str = "F"
    ) -> npt.NDArray[float]:
        """Return a 1-D array of Observable values.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
            order:  Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(self._collect("value", *obsid, err=err), order=order)

    def get_weighted_values(self, *obsid: str | int, err: float | None = None) -> tuple:
        """Return the weighted values of observables.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return self._collect("weighted_value", *obsid, err=err)

    def get_flat_weighted_values(
        self, *obsid: str | int, err: float | None = None, order: str = "F"
    ) -> npt.NDArray[float]:
        """Return a 1-D array of Observable weighted values.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
            order:  Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(self._collect("weighted_value", *obsid, err=err), order=order)

    def get_deviations(self, *obsid: str | int, err: float | None = None) -> tuple:
        """Return the deviations from target values.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return self._collect("deviation", *obsid, err=err)

    def get_flat_deviations(
        self, *obsid: str | int, err: float | None = None, order: str = "F"
    ) -> npt.NDArray[float]:
        """Return a 1-D array of deviations from target values.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
            order:  Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(self._collect("deviation", *obsid, err=err), order=order)

    def get_weighted_deviations(
        self, *obsid: str | int, err: float | None = None
    ) -> tuple:
        """Return the weighted deviations from target values.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return self._collect("weighted_deviation", *obsid, err=err)

    def get_flat_weighted_deviations(
        self, *obsid: str | int, err: float | None = None, order: str = "F"
    ) -> npt.NDArray[float]:
        """Return a 1-D array of weighted deviations from target values.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
            order:  Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(
            self._collect("weighted_deviation", *obsid, err=err), order=order
        )

    def get_weights(self, *obsid: str | int) -> tuple:
        """Return the weights of observables.

        Args:
            *obsid: name or index of selected observables (Default all)
        """
        return self._collect("weight", *obsid)

    def get_flat_weights(
        self, *obsid: str | int, order: str = "F"
    ) -> npt.NDArray[float]:
        """Return a 1-D array of Observable weights.

        Args:
            *obsid: name or index of selected observables (Default all)
            order:  Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(self._collect("weight", *obsid), order=order)

    def get_residuals(self, *obsid: str | int, err: float | None = None) -> tuple:
        """Return the residuals of observables.

        Args:
             *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return self._collect("residual", *obsid, err=err)

    def get_sum_residuals(self, *obsid: str | int, err: float | None = None) -> float:
        """Return the sum of residual values.

        Args:
            *obsid: name or index of selected observables (Default all)
            err:    Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return sum(np.sum(res) for res in self._collect("residual", *obsid, err=err))

    shapes = property(get_shapes, doc="Shapes of all values")
    flat_shape = property(get_flat_shape, doc="Shape of the flattened values")
    values = property(get_values, doc="values of all observables")
    flat_values = property(get_flat_values, doc="1-D array of Observable values")
    weighted_values = property(
        get_weighted_values, doc="Weighted values of all observables"
    )
    flat_weighted_values = property(
        get_flat_weighted_values, doc="1-D array of Observable weigthed values"
    )
    deviations = property(get_deviations, doc="Deviations from target values")
    flat_deviations = property(
        get_flat_deviations, doc="1-D array of deviations from target value"
    )
    weighted_deviations = property(
        get_weighted_deviations, doc="Weighted deviations from target values"
    )
    flat_weighted_deviations = property(
        get_flat_weighted_deviations,
        doc="1-D array of weighted deviations from target values",
    )
    flat_weights = property(get_flat_weights, doc="1-D array of Observable weights")
    residuals = property(get_residuals, doc="Residuals of all observable")
    sum_residuals = property(get_sum_residuals, doc="Sum of all residual values")
