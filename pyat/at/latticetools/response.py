import numpy as np
import copy
import abc
import multiprocessing
from functools import partial
from abc import ABC
import sys
from typing import Optional
if sys.version_info.minor < 9:
    from typing import Sequence
else:
    from collections.abc import Sequence
from ..lattice import AtError, Lattice, Refpts, Orbit
from .observables import ObservableList, OrbitObservable, RingObservable
from .observables import TrajectoryObservable
from .variables import ElementVariable, VariableList

_globring = None
_globobs = None


class SvdResponse(ABC):
    """SVD solver for corrections based on response matrices"""
    def __init__(self):
        self.weighted_response = None
        self.correction = None
        self.obsweights = None
        self.varweights = None
        self.v = None
        self.singular_values = None
        self.uh = None

    @abc.abstractmethod
    def build(self):
        ...

    def solve(self):
        """Computes the singular values of the response matrix"""
        if self.weighted_response is None:
            raise AtError("No response matrix: run build() first")
        u, s, vh = np.linalg.svd(self.weighted_response, full_matrices=False)
        self.v = vh.T * (1/s) * self.varweights.reshape(-1, 1)
        self.uh = u.T / self.obsweights
        self.singular_values = s

    def svd_solution(self, observed, nvals: Optional[int] = None):
        """"""
        if self.singular_values is None:
            self.solve()
        if nvals is None:
            nvals = len(self.singular_values)
        corr = self.v[:, :nvals] @ self.uh[:nvals, :] @ observed
        return corr

    def check_norm(self):
        """Displays the norm of the rows and columns of the weighted
        response matrix
        """
        if self.weighted_response is None:
            raise AtError("No response matrix: run build() first")
        obs = np.linalg.norm(self.weighted_response, axis=1)
        var = np.linalg.norm(self.weighted_response, axis=0)
        print("Observables: {}".format(np.amax(obs)/np.amin(obs)))
        print("Variables: {}".format(np.amax(var)/np.amin(var)))
        return obs, var

    def get_response(self, weighted: bool = False):
        """Return the response matrix

        Args:
            weighted:   If :py:obj:`True`, returns the weighted response
              matrix

        Returns:
            response:   Response matrix
        """
        if self.weighted_response is None:
            raise AtError("No response matrix: run build() first")
        if weighted:
            return self.weighted_response
        else:
            weights = self.obsweights.reshape(-1, 1) / self.varweights
            return self.weighted_response * weights

    def get_correction(self, nvals: Optional[int] = None):
        """Returns the correction matrix (pseudo-inverse of the response
        matrix)

        Args:
            nvals:  Desired number of singular values. If :py:obj:`None`, use
              all singular values

        Returns:
            corr:   Correction matrix
        """
        if self.singular_values is None:
            raise AtError("No correction computed yet: run solve() first")
        if nvals is None:
            nvals = len(self.singular_values)
        return self.v[:, :nvals] @ self.uh[:nvals, :]


class ResponseMatrix(SvdResponse):
    """Base class for response matrices"""
    def __init__(self, ring: Lattice,
                 variables: VariableList,
                 observables: ObservableList, *,
                 r_in: Orbit = None):
        """

        Args:
            ring:           Design lattice, used to compute the response
            variables:      List of :py:class:`.Variable`\ s
            observables:    List of :py:class:`.Observable`\s
            r_in:
        """
        self.ring = ring
        self.variables = variables
        self.observables = observables
        self.r_in = r_in
        super().__init__()

    def __iadd__(self, other: "ResponseMatrix"):
        if not isinstance(other, type(self)):
            mess = "Cannot add a {} to an {}}"
            raise TypeError(mess.format(type(other), type(self)))
        self.variables += other.variables
        self.observables += other.observables
        self.response = None
        self.correction = None
        return self

    def __add__(self, other: "ResponseMatrix"):
        nresp = copy.copy(self)
        nresp += other
        return nresp

    def correct(self, ring: Lattice, nvals: int = None,
                niter: int = 1, apply: bool = False):
        """Compute and optionally apply the correction

        Args:
            ring:       Lattice to be corrected
            nvals:      Desired number of singular values. If :py:obj:`None`,
              use all singular values
            apply:      If :py:obj:`True`, apply the correction to *ring*
            niter:      Number of iterations. For more than one iteration,
              *apply* must be :py:obj:`True`

        Returns:

        """
        for _ in range(niter):
            self.observables.evaluate(ring, r_in=self.r_in)
            corr = self.svd_solution(self.observables.flat_values, nvals=nvals)
            if apply:
                self.variables.increment(ring, corr)
        return corr

    @staticmethod
    def _resp_one(ring, observables, variable):
        if ring is None:
            ring = _globring
        if observables is None:
            observables = _globobs
        variable.step_up(ring)
        observables.evaluate(ring)
        op = observables.flat_weighted_values
        variable.step_down(ring)
        observables.evaluate(ring)
        om = observables.flat_weighted_values
        variable.set_initial(ring)
        return 0.5 * (op - om)

    def build(self, use_mp: bool = False, pool_size: Optional[int] = None,
              start_method: Optional[str] = None) -> None:
        """Build the response matrix

        Args:
            use_mp:             Use multiprocessing
            pool_size:          number of processes. If None,
              :pycode:`min(len(self.variables, nproc)` is used
            start_method:       python multiprocessing start method.
              :py:obj:`None` uses the python default that is considered safe.
              Available values: ``'fork'``, ``'spawn'``, ``'forkserver'``.
              Default for linux is ``'fork'``, default for macOS and  Windows
              is ``'spawn'``. ``'fork'`` may be used on macOS to speed up the
              calculation or to solve Runtime Errors, however it is considered
              unsafe.
        """
        boolrefs = self.ring.get_bool_index(None)
        for var in self.variables:
            boolrefs |= self.ring.get_bool_index(var.refpts)
            var.get(self.ring, initial=True)

        ring = self.ring.replace(boolrefs)

        if use_mp:
            ctx = multiprocessing.get_context(start_method)
            if pool_size is None:
                pool_size = min(len(self.variables),
                                multiprocessing.cpu_count())
            if ctx.get_start_method() == 'fork':
                global _globring
                global _globobs
                _globring = ring
                _globobs = self.observables
                args = [(None, None, var) for var in self.variables]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
                _globring = None
                _globobs = None
            else:
                args = [(ring, self.observables, var)
                        for var in self.variables]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
        else:
            results = [self._resp_one(ring, self.observables, var)
                       for var in self.variables]
        self.obsweights = self.observables.flat_weights
        self.varweights = self.variables.delta
        self.weighted_response = np.stack(results, axis=-1)

    def exclude(self, obsname: str, excluded: Refpts) -> None:
        """Exclude items from :py:class:`.Observable`\ s

        Args:
            obsname:    :py:class:`.Observable` name.
            excluded:   location of elements to excluse

        Example:
            >>> resp = OrbitResponseMatrix(ring, Monitor, Corrector)
            >>> obslist.exclude('Orbit[0]', 'BPM_02')

            Create an :py:class:`OrbitResponseMatrix` from
            :py:class:`.Corrector` elements to :py:class:`.Monitor` elements,
            and exclude the monitor with name "BPM_02"
        """
        self.observables.exclude(obsname, excluded)


class OrbitResponseMatrix(ResponseMatrix):
    """Orbit response matrix"""
    def __init__(self, ring: Lattice,  plane: int,
                 bpmrefs: Refpts,
                 steerrefs: Refpts, *,
                 cavrefs: Refpts = None,
                 bpmweight: float = 1.0,
                 steerdelta: float = 0.0001,
                 cavdelta: float = 100.0,
                 steersum: bool = False):
        """

        Args:
            ring:
            plane:
            bpmrefs:
            steerrefs:
            cavrefs:
            bpmweight:
            steerdelta:
            cavdelta:
            steersum:
        """
        def steerer(idb):
            return ElementVariable(idb, 'KickAngle', index=plane,
                                   delta=steerdelta)
        # Observables
        bpms = OrbitObservable(bpmrefs, axis=2*plane, weight=bpmweight)
        observables = ObservableList(ring, [bpms])
        if steersum:
            observables.append(RingObservable(steerrefs, 'KickAngle',
                                              index=plane, statfun=np.sum))
        # Variables
        steerers = (steerer(idx) for idx in ring.get_uint32_index(steerrefs))
        variables = VariableList(steerers)
        if cavrefs is not None:
            active = (el.longt_motion for el in ring.select(cavrefs))
            if not all(active):
                raise ValueError("Cavities are not active")
            variables.append(ElementVariable(cavrefs, 'Frequency',
                                             delta=cavdelta))

        super().__init__(ring, variables, observables)


class TrajectoryResponseMatrix(ResponseMatrix):
    """Trajectory response matrix"""
    def __init__(self, ring: Lattice, plane: int,
                 bpmrefs: Refpts,
                 steerrefs: Refpts, *,
                 r_in: Orbit = None,
                 bpmweight: float = 1.0,
                 steerdelta: float = 0.0001,
                 steersum: bool = False):
        """

        Args:
            ring:
            plane:
            bpmrefs:
            steerrefs:
            r_in:
            bpmweight:
            steerdelta:
            steersum:
        """
        def steerer(idb):
            return ElementVariable(idb, 'KickAngle', index=plane,
                                   delta=steerdelta)

        # Observables
        bpms = TrajectoryObservable(bpmrefs, axis=2*plane, weight=bpmweight)
        observables = ObservableList(ring, [bpms])
        if steersum:
            observables.append(RingObservable(steerrefs, 'KickAngle',
                                              index=plane, statfun=np.sum))
        # Variables
        steerers = (steerer(idx) for idx in ring.get_uint32_index(steerrefs))
        variables = VariableList(steerers)

        super().__init__(ring, variables, observables, r_in=r_in)
