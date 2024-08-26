r"""

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

from collections.abc import Iterable
from functools import reduce
from typing import Optional, Callable

import numpy as np

# noinspection PyProtectedMember
from .observables import Observable, ElementObservable, Need
from ..lattice import AtError, frequency_control
from ..lattice import Lattice, Orbit, Refpts, All
from ..physics import linopt6
from ..tracking import internal_lpass


def _flatten(vals, order="F"):
    return np.concatenate([np.reshape(v, -1, order=order) for v in vals])


class ObservableList(list):
    """Handles a list of Observables to be evaluated together

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

            >>> obslist.get_flat_values('tune', 'emittances[h]')
            array([3.815630e-01, 8.543764e-01, 1.090607e-04, 1.320391e-10])

            Get a flattened array of tunes and horizontal emittance
        """
        self.orbitrefs = None
        self.opticsrefs = None
        self.passrefs = None
        self.matrixrefs = None
        self.needs = None
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
        for obs in self:
            needs |= obs.needs
            if isinstance(obs, ElementObservable):
                needs.add(Need.RING)

        if (needs & self.needs_ring) and ring is None:
            raise ValueError("At least one Observable needs a ring argument")

        self.needs = needs
        if ring is None:
            # Initialise each observable
            for obs in self:
                obs._setup(ring)
        else:
            # Initialise each observable and make a summary all refpoints
            noref = ring.get_bool_index(None)
            orbitrefs = opticsrefs = passrefs = matrixrefs = noref
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
            self.orbitrefs = orbitrefs
            self.opticsrefs = opticsrefs
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
        if not isinstance(obs, Observable):
            raise TypeError(f"Cannot append a {type(obs)} to an ObservableList")
        self.needs = None
        super().append(obs)

    def extend(self, obsiter: Iterable[Observable]):
        self.needs = None
        super().extend(obsiter)

    def insert(self, index: int, obs: Observable):
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
        ring: Optional[Lattice] = None,
        *,
        dp: Optional[float] = None,
        dct: Optional[float] = None,
        df: Optional[float] = None,
        initial: bool = False,
        **kwargs,
    ):
        r"""Compute all the :py:class:`Observable` values

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
            """Evaluate a single observable"""

            def check_error(data, refpts):
                return data if isinstance(data, AtError) else data[refpts]

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
            obs.evaluate(*data, initial=initial)

        @frequency_control
        def ringeval(
            ring,
            dp: Optional[float] = None,
            dct: Optional[float] = None,
            df: Optional[float] = None,
        ):
            """Optics computations"""
            trajs = orbits = rgdata = eldata = emdata = mxdata = geodata = None
            o0 = None
            needs = self.needs

            if Need.TRAJECTORY in needs:
                # Trajectory computation
                r_in = kwargs.get("r_in", self.r_in)
                if r_in is None:
                    r_in = np.zeros(6)
                r_out = internal_lpass(ring, r_in.copy(), 1, refpts=self.passrefs)
                trajs = r_out[:, 0, :, 0].T

            if needs & self.needs_orbit:
                # Closed orbit computation
                orbit0 = kwargs.get("orbit", self.orbit)
                try:
                    o0, orbits = ring.find_orbit(
                        refpts=self.orbitrefs, dp=dp, dct=dct, df=df, orbit=orbit0
                    )
                except AtError as err:
                    orbits = mxdata = rgdata = eldata = emdata = err

            if Need.MATRIX in needs and o0 is not None:
                # Transfer matrix computation
                if ring.is_6d:
                    # noinspection PyUnboundLocalVariable
                    _, mxdata = ring.find_m66(
                        refpts=self.matrixrefs,
                        dp=dp,
                        dct=dct,
                        df=df,
                        orbit=o0,
                        keep_lattice=True,
                    )
                else:
                    # noinspection PyUnboundLocalVariable
                    _, mxdata = ring.find_m44(
                        refpts=self.matrixrefs,
                        dp=dp,
                        dct=dct,
                        df=df,
                        orbit=o0,
                        keep_lattice=True,
                    )

            if (needs & self.needs_optics) and o0 is not None:
                # Linear optics computation
                try:
                    _, rgdata, eldata = ring.get_optics(
                        refpts=self.opticsrefs,
                        dp=dp,
                        dct=dct,
                        df=df,
                        orbit=o0,
                        keep_lattice=True,
                        get_chrom=Need.CHROMATICITY in needs,
                        get_w=Need.W_FUNCTIONS in needs,
                        twiss_in=kwargs.get("twiss_in", self.twiss_in),
                        method=kwargs.get("method", self.method),
                    )
                except AtError as err:
                    rgdata = eldata = err

            if Need.EMITTANCE in needs and o0 is not None:
                # Emittance computation
                try:
                    emdata = ring.envelope_parameters(orbit=o0, keep_lattice=True)
                except Exception as err:
                    emdata = err

            if Need.GEOMETRY in needs:
                # Geometry computation
                geodata, _ = ring.get_geometry()

            return trajs, orbits, rgdata, eldata, emdata, mxdata, geodata

        if self.needs is None or initial:
            self._setup(ring)

        trajs, orbits, rgdata, eldata, emdata, mxdata, geodata = ringeval(
            ring, dp=dp, dct=dct, df=df
        )
        for ob in self:
            obseval(ring, ob)

    # noinspection PyProtectedMember
    def exclude(self, obsname: str, excluded: Refpts):
        # Set the excluded mask on the selected observable
        for obs in self:
            if obs.name == obsname:
                obs._excluded = excluded
        self.needs = None

    def _select(self, *obsnames: str):
        """Return an iterable over selected observables"""
        if obsnames:
            sel = set(obsnames)
            return (obs for obs in self if obs.name in sel)
        else:
            return self

    def _substitute(self, attrname: str, *obsnames: str, err: Optional[float] = None):
        for obs in self._select(*obsnames):
            try:
                res = getattr(obs, attrname)
            except AtError:
                if err is None:
                    raise
                else:
                    # noinspection PyProtectedMember
                    shp = obs._shape
                    res = np.broadcast_to(err, () if shp is None else shp)
            yield res

    def get_shapes(self, *obsnames: str) -> list:
        """Return the shapes of all values

        Args:
            *obsnames: names of selected observables (Default all)
        """
        return list(self._substitute("_shape", *obsnames))

    def get_flat_shape(self, *obsnames: str):
        """Shape of the flattened values

        Args:
            *obsnames: names of selected observables (Default all)
        """
        vals = (
            reduce(lambda x, y: x * y, shp, 1)
            for shp in self._substitute("_shape", *obsnames)
        )
        return (sum(vals),)

    def get_values(self, *obsnames: str, err: Optional[float] = None) -> list:
        """Return the values of observables

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return list(self._substitute("value", *obsnames, err=err))

    def get_flat_values(
        self, *obsnames: str, err: Optional[float] = None, order: str = "F"
    ) -> np.ndarray:
        """Return a 1-D array of Observable values

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(self._substitute("value", *obsnames, err=err), order=order)

    def get_weighted_values(self, *obsnames: str, err: Optional[float] = None) -> list:
        """Return the weighted values of observables

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return list(self._substitute("weighted_value", *obsnames, err=err))

    def get_flat_weighted_values(
        self, *obsnames: str, err: Optional[float] = None, order: str = "F"
    ) -> np.ndarray:
        """Return a 1-D array of Observable weighted values

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(
            self._substitute("weighted_value", *obsnames, err=err), order=order
        )

    def get_deviations(self, *obsnames: str, err: Optional[float] = None) -> list:
        """Return the deviations from target values

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return list(self._substitute("deviation", *obsnames, err=err))

    def get_flat_deviations(
        self, *obsnames: str, err: Optional[float] = None, order: str = "F"
    ) -> np.ndarray:
        """Return a 1-D array of deviations from target values

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(self._substitute("deviation", *obsnames, err=err), order=order)

    def get_weighted_deviations(
        self, *obsnames: str, err: Optional[float] = None
    ) -> list:
        """Return the weighted deviations from target values

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return list(self._substitute("weighted_deviation", *obsnames, err=err))

    def get_flat_weighted_deviations(
        self, *obsnames: str, err: Optional[float] = None, order: str = "F"
    ) -> np.ndarray:
        """Return a 1-D array of weighted deviations from target values

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(
            self._substitute("weighted_deviation", *obsnames, err=err), order=order
        )

    def get_weights(self, *obsnames: str) -> list:
        """Return the weights of observables

        Args:
            *obsnames: names of selected observables (Default all)
        """
        return list(self._substitute("weight", *obsnames))

    def get_flat_weights(self, *obsnames: str, order: str = "F") -> np.ndarray:
        """Return a 1-D array of Observable weights

        Args:
            *obsnames: names of selected observables (Default all)
            order:      Ordering for reshaping. See :py:func:`~numpy.reshape`
        """
        return _flatten(self._substitute("weight", *obsnames), order=order)

    def get_residuals(self, *obsnames: str, err: Optional[float] = None) -> list:
        """Return the residuals of observables

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return list(self._substitute("residual", *obsnames, err=err))

    def get_sum_residuals(self, *obsnames: str, err: Optional[float] = None) -> float:
        """Return the sum of residual values

        Args:
            *obsnames: names of selected observables (Default all)
            err: Default observable value to be used when the evaluation failed. By
              default, an Exception is raised.
        """
        return sum(
            np.sum(res) for res in self._substitute("residual", *obsnames, err=err)
        )

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
