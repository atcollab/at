"""
Definition of :py:class:`.ResponseMatrix` objects
"""
from __future__ import annotations
import numpy as np
import copy
import abc
import multiprocessing
from functools import partial
from abc import ABC
from typing import Optional
from ..lattice import AtError, Lattice, Refpts, Orbit, AxisDef, plane_
from .observables import ObservableList, OrbitObservable, LatticeObservable
from .observables import TrajectoryObservable
from .variables import Variable, ElementVariable, VariableList

_globring = None
_globobs = None


def _resp_one(ring: Lattice, observables: ObservableList, variable: Variable):
    variable.step_up(ring)
    observables.evaluate(ring)
    op = observables.flat_weighted_values
    variable.step_down(ring)
    observables.evaluate(ring)
    om = observables.flat_weighted_values
    variable.set_initial(ring)
    return 0.5 * (op - om)


def _resp_one_fork(variable: Variable):
    # noinspection PyTypeChecker
    return _resp_one(_globring, _globobs, variable)


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
        self._obsmask = None
        self._varmask = None

    @abc.abstractmethod
    def build(self):
        nobs, nvar = self.weighted_response.shape
        self._obsmask = np.ones(nobs, dtype=bool)
        self._varmask = np.ones(nvar, dtype=bool)

    def solve(self):
        """Compute the singular values of the response matrix"""
        resp = self.weighted_response
        if resp is None:
            raise AtError("No response matrix: run build() first")
        selected = np.ix_(self._obsmask, self._varmask)
        u, s, vh = np.linalg.svd(resp[selected], full_matrices=False)
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
        """Display the norm of the rows and columns of the weighted
        response matrix

        Adjusting the variables and observable weights to equalize the norms
        of rows and columns is important.

        Returns:
            obs_norms:      Norms of observables (rows)
            var_norms:      Norms of Variables (columns)
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
        """Return the correction matrix (pseudo-inverse of the response
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
    r"""Base class for response matrices

    It is defined by any arbitrary set of :py:class:`.Variable`\ s and
    :py:class:`.Observable`\s

    Addition is defined on :py:class:`ResponseMatrix` objects as the addition
    of their :py:class:`.Variable`\ s and :py:class:`.Observable`\s to produce
    combined responses"""
    def __init__(self, variables: VariableList, observables: ObservableList, *,
                 r_in: Orbit = None):
        r"""

        Args:
            variables:      List of :py:class:`.Variable`\ s
            observables:    List of :py:class:`.Observable`\s
            r_in:
        """
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

    def correct(self, ring: Lattice,
                nvals: int = None,
                niter: int = 1, apply: bool = False):
        """Compute and optionally apply the correction

        Args:
            ring:       Lattice description. The response matrix observables
              will be evaluated for *ring* and the deviation from target will
              be   corrected
            nvals:      Desired number of singular values. If :py:obj:`None`,
              use all singular values
            apply:      If :py:obj:`True`, apply the correction to *ring*
            niter:      Number of iterations. For more than one iteration,
              *apply* must be :py:obj:`True`

        Returns:
            correction: Vector of correction values
        """
        if niter > 1 and not apply:
            raise ValueError("needs: appli is True")
        obs = self.observables
        if apply:
            self.variables.get(ring)
        sumcorr = np.array([0.0])
        for _ in range(niter):
            obs.evaluate(ring, r_in=self.r_in)
            corr = self.svd_solution(-obs.flat_deviations, nvals=nvals)
            sumcorr = sumcorr + corr    # non-broadcastable sumcorr
            if apply:
                self.variables.increment(ring, corr)
        return sumcorr

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
        ring = self.observables.ring
        boolrefs = ring.get_bool_index(None)
        for var in self.variables:
            boolrefs |= ring.get_bool_index(var.refpts)
            var.get(ring, initial=True)

        ring = ring.replace(boolrefs)
        self.observables.evaluate(ring)
        self.obsweights = self.observables.flat_weights
        self.varweights = self.variables.deltas

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
                with ctx.Pool(pool_size) as pool:
                    results = pool.map(_resp_one_fork, self.variables)
                _globring = None
                _globobs = None
            else:
                _resp_one_spawn = partial(_resp_one, ring, self.observables)
                with ctx.Pool(pool_size) as pool:
                    results = pool.map(_resp_one_spawn, self.variables)
        else:
            results = [_resp_one(ring, self.observables, var)
                       for var in self.variables]
        self.weighted_response = np.stack(results, axis=-1)
        super().build()

    def exclude_obs(self, obsname: str, excluded: Refpts) -> None:
        # noinspection PyUnresolvedReferences
        r"""Exclude items from :py:class:`.Observable`\ s

        Args:
            obsname:    :py:class:`.Observable` name.
            excluded:   location of elements to excluse

        Example:
            >>> resp = OrbitResponseMatrix(ring, 'h', Monitor, Corrector)
            >>> resp.exclude_obs('h_orbit', 'BPM_02')

            Create an horizontal :py:class:`OrbitResponseMatrix` from
            :py:class:`.Corrector` elements to :py:class:`.Monitor` elements,
            and exclude the monitor with name "BPM_02"
        """
        self.observables.exclude(obsname, excluded)
        # Force a full rebuild
        self.singular_values = None
        self.weighted_response = None


class OrbitResponseMatrix(ResponseMatrix):
    r"""Orbit response matrix

    An :py:class:`OrbitResponseMatrix` applies to a single plane, horizontal or
    vertical. A combined response matrix is obtained by adding horizontal and
    vertical matrices.

    Variables are a set of steerers and optionally the RF frequency. Steerer
    variables are named ``hxxxx`` or ``vxxxx`` where xxxx is the index in the
    lattice. The RF frequency variable is named ``RF frequency``.

    Observables are the closed orbit position at selected points, named
    ``h_orbit`` for the horizontal plane or ``v_orbit`` for vertical plane,
    and optionally the sum of steerer angles named ``sum(h_kicks)`` or
    ``sum(v_kicks)``
    """
    def __init__(self, ring: Lattice,  plane: AxisDef,
                 bpmrefs: Refpts,
                 steerrefs: Refpts, *,
                 cavrefs: Refpts = None,
                 bpmweight: float = 1.0,
                 bpmtarget=0.0,
                 steerdelta: float = 0.0001,
                 cavdelta: float = 100.0,
                 steersum: bool = False):
        """

        Args:
            ring:       Design lattice, used to compute the response
            plane:      One out of {0, 'x', 'h', 'H'} for horizontal orbit, or
              one of {1, 'y', 'v', 'V'} for vertical orbit
            bpmrefs:    Location of closed orbit observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            steerrefs:  Location of orbit steerers. Their *KickAngle* attribute
              is used.
            cavrefs:    Location of RF cavities. Their *Frequency* attribute
              is used. If :py:obj:`None`, no cavity is included in the response.
              Cavities must be active.
            bpmweight:  Weight on position readings
            bpmtarget:  Target position.
            steerdelta: Step on steerers for matrix computation [rad]. This is
              also the steerer weight
            cavdelta:   Step on RF frequency for matrix computation [Hz]. This
              is also the cavity weight
            steersum:   If :py:obj:`True`, the sum of steerers is added to the
              Observables
        """
        def steerer(ik):
            name = f"{plcode}{ik:04}"
            return ElementVariable(ik, 'KickAngle', index=pl, name=name,
                                   delta=steerdelta)

        pl = plane_(plane, 'index')
        plcode = plane_(plane, 'code')
        # Observables
        nm = f"{plcode}_orbit"
        bpms = OrbitObservable(bpmrefs, axis=2*pl, name=nm, target=bpmtarget,
                               weight=bpmweight)
        observables = ObservableList(ring, [bpms])
        if steersum:
            nm = f"{plcode}_kicks"
            observables.append(LatticeObservable(steerrefs, 'KickAngle',
                                                 name=nm,
                                                 target=0.0, index=pl,
                                                 statfun=np.sum))
        # Variables
        steerers = (steerer(idx) for idx in ring.get_uint32_index(steerrefs))
        variables = VariableList(steerers)
        if cavrefs is not None:
            active = (el.longt_motion for el in ring.select(cavrefs))
            if not all(active):
                raise ValueError("Cavities are not active")
            variables.append(ElementVariable(cavrefs, 'Frequency',
                                             name="RF frequency",
                                             delta=cavdelta))

        super().__init__(variables, observables)


class TrajectoryResponseMatrix(ResponseMatrix):
    """Trajectory response matrix

    Variables are a set of steerers,

    Observables are the trajectory position at selected points, named
    ``h_positions`` for the horizontal plane or ``v_positions`` for vertical
    plane.

    An :py:class:`TrajectoryResponseMatrix` applies to a single plane,
    horizontal or vertical. A combined response matrix is obtained by adding
    horizontal and vertical matrices"""

    def __init__(self, ring: Lattice, plane: AxisDef,
                 bpmrefs: Refpts,
                 steerrefs: Refpts, *,
                 r_in: Orbit = None,
                 bpmweight: float = 1.0,
                 bpmtarget=0.0,
                 steerdelta: float = 0.0001,
                 steersum: bool = False):
        """

        Args:
            ring:       Design lattice, used to compute the response
            plane:      One out of {0, 'x', 'h', 'H'} for horizontal orbit, or
              one of {1, 'y', 'v', 'V'} for vertical orbit
            bpmrefs:    Location of closed orbit observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            steerrefs:  Location of orbit steerers. Their *KickAngle* attribute
              is used.
            r_in:       (6,) vector of initial coordinates of the trajectory
            bpmweight:  Weight on position readings
            bpmtarget:  Target position
            steerdelta: Step on steerers for matrix computation [rad]. This is
              also the steerer weight
        """
        def steerer(ik):
            name = f"{plcode}{ik:04}"
            return ElementVariable(ik, 'KickAngle', index=pl, name=name,
                                   delta=steerdelta)

        pl = plane_(plane, 'index')
        plcode = plane_(plane, 'code')
        # Observables
        nm = f"{plcode}_positions"
        bpms = TrajectoryObservable(bpmrefs, axis=2*pl, name=nm,
                                    target=bpmtarget,
                                    weight=bpmweight)
        observables = ObservableList(ring, [bpms])
        if steersum:
            nm = f"{plcode}_kicks"
            observables.append(LatticeObservable(steerrefs, 'KickAngle',
                                                 name=nm,
                                                 target=0.0, index=plane,
                                                 statfun=np.sum))
        # Variables
        steerers = (steerer(idx) for idx in ring.get_uint32_index(steerrefs))
        variables = VariableList(steerers)

        super().__init__(variables, observables, r_in=r_in)
