# noinspection PyUnresolvedReferences
r"""Definition of :py:class:`.ResponseMatrix` objects.

A :py:class:`ResponseMatrix` object defines a general-purpose response matrix, based
on a :py:class:`.VariableList` of attributes which will be independently varied, and an
:py:class:`.ObservableList` of attributes which will be recorded for each
variable step.

:py:class:`ResponseMatrix` objects can be combined with the "+" operator to define
combined responses. This concatenates the variables and the observables.

This module also defines two commonly used response matrices:
:py:class:`OrbitResponseMatrix` for circular machines and
:py:class:`TrajectoryResponseMatrix` for beam lines. Other matrices can be easily
defined by providing the desired Observables and Variables to the
:py:class:`ResponseMatrix` base class.

Instantiation
^^^^^^^^^^^^^

The simplest orbit response matrix can be instantiated with:

>>> resp_v = OrbitResponseMatrix(ring, "v")

By default, the observables are all the :py:class:`.Monitor` elements, and the
variables are all the elements having a *KickAngle* attribute. This is equivalent to:

>>> resp_v = OrbitResponseMatrix(
...     ring, "v", bpmrefs=at.Monitor, steerrefs=at.checkattr("KickAngle")
... )

If correction is desired, the variable elements must have the *KickAngle* attribute
used for correction. It's available for all magnets, though not present by default
except in :py:class:`.Corrector` magnets. For other magnets, the attribute should be
explicitly created.

There are options to include the RF frequency in the variable list, and the sum of
correction angles in the list of observables:

>>> resp_h = OrbitResponseMatrix(ring, "h", cavrefs=at.RFCavity, steersum=True)

A combined horizontal+vertical response matrix is obtained with:

>>> resp_hv = resp_h + resp_v

Matrix Building
^^^^^^^^^^^^^^^

The response matrix may be built by two methods:

1. :py:meth:`~ResponseMatrix.build` computes the matrix using tracking.
2. :py:meth:`~ResponseMatrix.load` loads data from a file containing previously
   saved values or experimentally measured values

Normalisation
^^^^^^^^^^^^^

To be correctly inverted, the response matrix must be correctly normalised: the norms
of its columns must be of the same order of magnitude, and similarly for the rows.
This is done by adjusting the weights :math:`w_v` for the variables :math:`\mathbf{V}`
and :math:`w_o` for the observables :math:`\mathbf{O}`.
With :math:`\mathbf{R}` the response matrix:

.. math::

   \mathbf{O} = \mathbf{R} . \mathbf{V}

The weighted response matrix :math:`\mathbf{R}_w` is:

.. math::

   \frac{\mathbf{O}}{w_o} = \mathbf{R}_w . \frac{\mathbf{V}}{w_v}

The :math:`\mathbf{R}_w` is dimensionless and should be normalised. This can be checked
using:

* :py:meth:`~ResponseMatrix.check_norm` which prints the ratio of the maximum / minimum
  norms for variables and observables. These should be less than 10.
* :py:meth:`~.ResponseMatrix.plot_norm`

Both natural and weighted response matrices can be retrieved with the
:py:attr:`~ResponseMatrix.response` and :py:attr:`~ResponseMatrix.weighted_response`
properties.

Inversion
^^^^^^^^^

The :py:meth:`~ResponseMatrix.solve` method computes the singular values of the
weighted response matrix.

After solving, orbit correction is available, for instance with

* :py:meth:`~ResponseMatrix.get_correction` which returns the correction matrix,
* :py:meth:`~ResponseMatrix.correct` which computes and optionally applies a correction
  for the provided :py:class:`.Lattice`.

.. include:: ../notebooks/ATPrimer.rst


"""

from __future__ import annotations

__all__ = [
    "sequence_split",
    "ResponseMatrix",
    "OrbitResponseMatrix",
    "TrajectoryResponseMatrix",
]

import os
import copy
import multiprocessing
import concurrent.futures
import abc
from collections.abc import Sequence, Generator
from itertools import chain
from functools import partial
import math

import numpy as np

from .observables import TrajectoryObservable, OrbitObservable, LatticeObservable
from .observables import LocalOpticsObservable, GlobalOpticsObservable
from .observablelist import ObservableList
from ..lattice import AtError, Lattice, Refpts, AxisDef, plane_
from ..lattice import Monitor, checkattr
from ..lattice.lattice_variables import RefptsVariable
from ..lattice.variables import VariableList

_orbit_correctors = checkattr("KickAngle")

_globring: Lattice | None = None
_globobs: ObservableList | None = None


def sequence_split(seq: Sequence, nslices: int) -> Generator[Sequence, None, None]:
    """Split a sequence into multiple sub-sequences.

    The length of *seq* does not have to be a multiple of *nslices*.

    Args:
        seq: sequence to split
        nslices: number of sub-sequences

    Returns:
        subseqs: Iterator over sub-sequences
    """

    def _split(seqsizes):
        beg = 0
        for size in seqsizes:
            end = beg + size
            yield seq[beg:end]
            beg = end

    lna = len(seq)
    sz, rem = divmod(lna, nslices)
    lsubseqs = [sz] * nslices
    for k in range(rem):
        lsubseqs[k] += 1
    return _split(lsubseqs)


def _resp(
    ring: Lattice, observables: ObservableList, variables: VariableList, **kwargs
):
    def _resp_one(variable: RefptsVariable):
        """Single response"""
        variable.step_up(ring=ring)
        observables.evaluate(ring, **kwargs)
        op = observables.flat_values
        variable.step_down(ring=ring)
        observables.evaluate(ring, **kwargs)
        om = observables.flat_values
        variable.reset(ring=ring)
        return (op - om) / (2.0 * variable.delta)

    return [_resp_one(v) for v in variables]


def _resp_fork(variables: VariableList, **kwargs):
    """Response for fork parallel method."""
    return _resp(_globring, _globobs, variables, **kwargs)


class _SvdSolver(abc.ABC):
    """SVD solver for response matrices."""

    def __init__(self):
        self._response = None
        self.v = None
        self.singular_values = None
        self.uh = None
        self._obsmask = None
        self._varmask = None

    @abc.abstractmethod
    def build(self) -> None:
        """Build the response matrix."""
        nobs, nvar = self._response.shape
        self._obsmask = np.ones(nobs, dtype=bool)
        self._varmask = np.ones(nvar, dtype=bool)

    @property
    @abc.abstractmethod
    def varweights(self): ...

    @property
    @abc.abstractmethod
    def obsweights(self): ...

    def solve(self) -> None:
        """Compute the singular values of the response matrix."""
        resp = self.weighted_response
        if resp is None:
            raise AtError("No response matrix: run build() first")
        selected = np.ix_(self._obsmask, self._varmask)
        u, s, vh = np.linalg.svd(resp[selected], full_matrices=False)
        self.v = vh.T * (1 / s) * self.varweights.reshape(-1, 1)
        self.uh = u.T / self.obsweights
        self.singular_values = s

    def check_norm(self) -> tuple[float, float]:
        """Display the norm of the rows and columns of the weighted response matrix.

        Adjusting the variables and observable weights to equalize the norms
        of rows and columns is important.

        Returns:
            obs_norms:      Norms of observables (rows)
            var_norms:      Norms of Variables (columns)
        """
        resp = self.weighted_response
        if resp is None:
            raise AtError("No response matrix: run build() first")
        obs = np.linalg.norm(resp, axis=1)
        var = np.linalg.norm(resp, axis=0)
        print(f"max/min Observables: {np.amax(obs) / np.amin(obs)}")
        print(f"max/min Variables: {np.amax(var) / np.amin(var)}")
        return obs, var

    @property
    def response(self):
        """Response matrix."""
        return self._response

    @property
    def weighted_response(self):
        """Weighted response matrix."""
        resp = self._response
        if resp is None:
            return None
        else:
            return resp * (self.varweights / self.obsweights.reshape(-1, 1))

    def correction_matrix(self, nvals: int | None = None):
        """Return the correction matrix (pseudo-inverse of the response matrix).

        Args:
            nvals:  Desired number of singular values. If :py:obj:`None`, use
              all singular values

        Returns:
            cormat: Correction matrix
        """
        if self.singular_values is None:
            self.solve()
        if nvals is None:
            nvals = len(self.singular_values)
        return self.v[:, :nvals] @ self.uh[:nvals, :]

    def get_correction(self, observed, nvals: int | None = None):
        """Compute the correction of the given observation.

        Args:
            observed:   Observed error vector
            nvals:      Desired number of singular values. If :py:obj:`None`, use
              all singular values

        Returns:
            corr:       Correction vector
        """
        return -self.correction_matrix(nvals=nvals) @ observed

    def save(self, file) -> None:
        """Save a response matrix.

        Args:
            file:   file-like object, string, or :py:class:`pathlib.Path`: File to
              which the data is saved. If file is a file-object, it must be opened in
              binary mode. If file is a string or Path, a .npy extension will
              be appended to the filename if it does not already have one.
        """
        if self._response is None:
            raise AtError("No response matrix: run build() first")
        np.save(file, self._response)

    def load(self, file) -> None:
        """Load a response matrix.

        Args:
            file:   file-like object, string, or :py:class:`pathlib.Path`: the file to
              read. A file object must always be opened in binary mode.
        """
        self._response = np.load(file)
        nobs, nvar = self._response.shape
        self._obsmask = np.ones(nobs, dtype=bool)
        self._varmask = np.ones(nvar, dtype=bool)


class ResponseMatrix(_SvdSolver):
    r"""Base class for response matrices.

    It is defined by any arbitrary set of :py:class:`~.variables.VariableBase`\ s and
    :py:class:`.Observable`\s

    Addition is defined on :py:class:`ResponseMatrix` objects as the addition
    of their :py:class:`~.variables.VariableBase`\ s and :py:class:`.Observable`\s to
    produce combined responses.
    """

    def __init__(
        self,
        ring: Lattice,
        variables: VariableList,
        observables: ObservableList,
    ):
        r"""
        Args:
            ring:           Design lattice, used to compute the response
            variables:      List of :py:class:`~.variables.VariableBase`\ s
            observables:    List of :py:class:`.Observable`\s
        """
        # for efficiency of parallel computation, the variable's refpts must be integer
        for var in variables:
            var.refpts = ring.get_uint32_index(var.refpts)
        self.ring = ring
        self.variables = variables
        self.observables = observables
        self.buildargs = {}
        variables.get(ring=ring, initial=True)
        observables.evaluate(ring=ring, initial=True)
        super().__init__()

    def __iadd__(self, other: ResponseMatrix):
        if not isinstance(other, type(self)):
            mess = "Cannot add a {} to an {}}"
            raise TypeError(mess.format(type(other), type(self)))
        self.variables += other.variables
        self.observables += other.observables
        self._response = None
        return self

    def __add__(self, other: ResponseMatrix):
        nresp = copy.copy(self)
        nresp += other
        return nresp

    def __str__(self):
        no, nv = self.shape
        return f"{type(self).__name__}({no} observables, {nv} variables)"

    @property
    def shape(self):
        """Shape of the response matrix."""
        return len(self.observables.flat_values), len(self.variables)

    @property
    def varweights(self):
        """Variable weights."""
        return self.variables.deltas

    @property
    def obsweights(self):
        """Observable weights."""
        return self.observables.flat_weights

    def correct(
        self, ring: Lattice, nvals: int = None, niter: int = 1, apply: bool = False
    ):
        """Compute and optionally apply the correction.

        Args:
            ring:       Lattice description. The response matrix observables
              will be evaluated for *ring* and the deviation from target will
              be corrected
            nvals:      Desired number of singular values. If :py:obj:`None`,
              use all singular values
            apply:      If :py:obj:`True`, apply the correction to *ring*
            niter:      Number of iterations. For more than one iteration,
              *apply* must be :py:obj:`True`

        Returns:
            correction: Vector of correction values
        """
        if niter > 1 and not apply:
            raise ValueError("needs: apply is True")
        obs = self.observables
        if apply:
            self.variables.get(ring=ring)
        sumcorr = np.array([0.0])
        for _ in range(niter):
            obs.evaluate(ring, **self.buildargs)
            corr = self.get_correction(obs.flat_deviations, nvals=nvals)
            sumcorr = sumcorr + corr  # non-broadcastable sumcorr
            if apply:
                self.variables.increment(corr, ring=ring)
        return sumcorr

    def build(
        self,
        use_mp: bool = False,
        pool_size: int | None = None,
        start_method: str | None = None,
        **kwargs,
    ) -> None:
        """Build the response matrix.

        Args:
            use_mp:             Use multiprocessing
            pool_size:          number of processes. If None,
              :pycode:`min(len(self.variables, nproc)` is used
            start_method:       python multiprocessing start method.
              :py:obj:`None` uses the python default that is considered safe.
              Available values: ``'fork'``, ``'spawn'``, ``'forkserver'``.
              Default for linux is ``'fork'``, default for macOS and  Windows
              is ``'spawn'``. ``'fork'`` may be used on macOS to speed up the
              calculation, however it is considered unsafe.

        Keyword Args:
            dp (float):     Momentum deviation. Defaults to :py:obj:`None`
            dct (float):    Path lengthening. Defaults to :py:obj:`None`
            df (float):     Deviation from the nominal RF frequency.
              Defaults to :py:obj:`None`
            r_in (Orbit):   Initial trajectory, used for
              :py:class:`TrajectoryResponseMatrix`, Default: zeros(6)
        """
        self.buildargs = kwargs
        self.observables.evaluate(self.ring)

        if use_mp:
            global _globring
            global _globobs
            ctx = multiprocessing.get_context(start_method)
            if pool_size is None:
                pool_size = min(len(self.variables), os.cpu_count())
            obschunks = sequence_split(self.variables, pool_size)
            if ctx.get_start_method() == "fork":
                _globring = self.ring
                _globobs = self.observables
                _single_resp = partial(_resp_fork, **kwargs)

            else:
                _single_resp = partial(_resp, self.ring, self.observables, **kwargs)
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=pool_size,
                mp_context=ctx,
            ) as pool:
                results = list(chain(*pool.map(_single_resp, obschunks)))
            # with ctx.Pool(pool_size) as pool:
            #     results = pool.map(_single_resp, self.variables)
            _globring = None
            _globobs = None
        else:
            ring = self.ring
            boolrefs = ring.get_bool_index(None)
            for var in self.variables:
                boolrefs |= ring.get_bool_index(var.refpts)

            ring = ring.replace(boolrefs)
            results = _resp(ring.deepcopy(), self.observables, self.variables, **kwargs)

        self._response = np.stack(results, axis=-1)
        super().build()

    def exclude_obs(self, obsname: str, excluded: Refpts) -> None:
        # noinspection PyUnresolvedReferences
        r"""Exclude items from :py:class:`.Observable`\ s.

        Args:
            obsname:    :py:class:`.Observable` name.
            excluded:   location of elements to excluse

        Example:
            >>> resp = OrbitResponseMatrix(ring, "h", Monitor, Corrector)
            >>> resp.exclude_obs("x_orbit", "BPM_02")

            Create an horizontal :py:class:`OrbitResponseMatrix` from
            :py:class:`.Corrector` elements to :py:class:`.Monitor` elements,
            and exclude the monitor with name "BPM_02"
        """
        self.observables.exclude(obsname, excluded)
        # Force a full rebuild
        self.singular_values = None
        self._response = None


class OrbitResponseMatrix(ResponseMatrix):
    r"""Orbit response matrix.

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

    def __init__(
        self,
        ring: Lattice,
        plane: AxisDef,
        bpmrefs: Refpts = Monitor,
        steerrefs: Refpts = _orbit_correctors,
        *,
        cavrefs: Refpts = None,
        bpmweight: float = 1.0,
        bpmtarget: float | Sequence[float] = 0.0,
        steerdelta: float | Sequence[float] = 0.0001,
        cavdelta: float | None = None,
        steersum: bool = False,
        stsumweight: float | None = None,
    ):
        """
        Args:
            ring:       Design lattice, used to compute the response.
            plane:      One out of {0, 'x', 'h', 'H'} for horizontal orbit, or
              one of {1, 'y', 'v', 'V'} for vertical orbit.
            bpmrefs:    Location of closed orbit observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`".
              Default: all :py:class:`.Monitor` elements.
            steerrefs:  Location of orbit steerers. Their *KickAngle* attribute
              is used and must be present in the selected elements.
              Default: All Elements having a *KickAngle* attribute.
            cavrefs:    Location of RF cavities. Their *Frequency* attribute
              is used. If :py:obj:`None`, no cavity is included in the response.
              Cavities must be active. Cavity variables are appended to the steerer
              variables.
            bpmweight:  Weight of position readings. Must be broadcastable to the
              number of BPMs.
            bpmtarget:  Target orbit position. Must be broadcastable to the number of
              observation points.
            cavdelta:   Step on RF frequency for matrix computation [Hz]. This
              is also the cavity weight
            steerdelta: Step on steerers for matrix computation [rad]. This is
              also the steerer weight. Must be broadcastable to the number of steerers.
            cavdelta:   Step on RF frequency for matrix computation [Hz]. This
              is also the cavity weight
            steersum:   If :py:obj:`True`, the sum of steerers is appended to the
              Observables.
            stsumweight: Weight on steerer summation. Default 1.0.

        :ivar VariableList variables: matrix variables
        :ivar ObservableList observables: matrix observables

        """

        def steerer(ik, delta):
            name = f"{plcode}{ik:04}"
            return RefptsVariable(ik, "KickAngle", index=pl, name=name, delta=delta)

        def set_norm():
            bpm = LocalOpticsObservable(bpmrefs, "beta", plane=pl)
            sts = LocalOpticsObservable(ids, "beta", plane=pl)
            dsp = LocalOpticsObservable(bpmrefs, "dispersion", plane=2 * pl)
            tun = GlobalOpticsObservable("tune", plane=pl)
            obs = ObservableList([bpm, sts, dsp, tun])
            result = obs.evaluate(ring=ring)
            alpha = ring.disable_6d(copy=True).get_mcf(0)
            freq = ring.get_rf_frequency(cavpts=cavrefs)
            nr = np.outer(
                np.sqrt(result[0]) / bpmweight, np.sqrt(result[1]) * steerdelta
            )
            vv = np.mean(np.linalg.norm(nr, axis=0))
            vo = np.mean(np.linalg.norm(nr, axis=1))
            korb = 0.25 * math.sqrt(2.0) / math.sin(math.pi * result[3])
            cavd = vv * korb * alpha * freq / np.linalg.norm(result[2] / bpmweight)
            stsw = np.linalg.norm(deltas) / vo / korb
            return cavd, stsw

        pl = plane_(plane, "index")
        plcode = plane_(plane, "code")
        ids = ring.get_uint32_index(steerrefs)
        nbsteers = len(ids)
        deltas = np.broadcast_to(steerdelta, nbsteers)
        if steersum and stsumweight is None or cavrefs and cavdelta is None:
            cavd, stsw = set_norm()

        # Observables
        nm = f"{plcode}_orbit"
        bpms = OrbitObservable(
            bpmrefs, axis=2 * pl, name=nm, target=bpmtarget, weight=bpmweight
        )
        observables = ObservableList([bpms])
        if steersum:
            sumobs = LatticeObservable(
                steerrefs,
                "KickAngle",
                name=f"{plcode}_kicks",
                target=0.0,
                index=pl,
                weight=stsumweight if stsumweight else stsw,
                statfun=np.sum,
            )
            observables.append(sumobs)

        # Variables
        variables = VariableList(steerer(ik, delta) for ik, delta in zip(ids, deltas))
        if cavrefs is not None:
            active = (el.longt_motion for el in ring.select(cavrefs))
            if not all(active):
                raise ValueError("Cavities are not active")
            cavvar = RefptsVariable(
                cavrefs,
                "Frequency",
                name="RF frequency",
                delta=cavdelta if cavdelta else cavd,
            )
            variables.append(cavvar)

        self.nbsteers = nbsteers

        super().__init__(ring, variables, observables)

    @property
    def bpmweight(self):
        """Weight of position readings."""
        return self.observables[0].weight

    @bpmweight.setter
    def bpmweight(self, value):
        self.observables[0].weight = value

    @property
    def stsumweight(self):
        """Weight of steerer summation."""
        return self.observables[1].weight

    @stsumweight.setter
    def stsumweight(self, value):
        self.observables[1].weight = value

    @property
    def steerdelta(self):
        """Step and weight on steerers."""
        return self.variables[: self.nbsteers].deltas

    @steerdelta.setter
    def steerdelta(self, value):
        self.variables[: self.nbsteers].deltas = value

    @property
    def cavdelta(self):
        """Step and weight on RF frequency deviation."""
        return self.variables[self.nbsteers].delta

    @cavdelta.setter
    def cavdelta(self, value):
        self.variables[self.nbsteers].delta = value


class TrajectoryResponseMatrix(ResponseMatrix):
    """Trajectory response matrix.

    Variables are a set of steerers,

    Observables are the trajectory position at selected points, named
    ``h_positions`` for the horizontal plane or ``v_positions`` for vertical
    plane.

    A :py:class:`TrajectoryResponseMatrix` applies to a single plane,
    horizontal or vertical. A combined response matrix is obtained by adding
    horizontal and vertical matrices
    """

    def __init__(
        self,
        ring: Lattice,
        plane: AxisDef,
        bpmrefs: Refpts = Monitor,
        steerrefs: Refpts = _orbit_correctors,
        *,
        bpmweight: float = 1.0,
        bpmtarget: float = 0.0,
        steerdelta: float = 0.0001,
    ):
        """
        Args:
            ring:       Design lattice, used to compute the response
            plane:      One out of {0, 'x', 'h', 'H'} for horizontal orbit, or
              one of {1, 'y', 'v', 'V'} for vertical orbit
            bpmrefs:    Location of closed orbit observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`".
              Default: all :py:class:`.Monitor` elements.
            steerrefs:  Location of orbit steerers. Their *KickAngle* attribute
              is used and must be present in the selected elements.
              Default: All Elements having a *KickAngle* attribute.
            bpmweight:  Weight on position readings. Must be broadcastable to the
              number of BPMs
            bpmtarget:  Target position
            steerdelta: Step on steerers for matrix computation [rad]. This is
              also the steerer weight. Must be broadcastable to the number of steerers.
        """

        def steerer(ik):
            name = f"{plcode}{ik:04}"
            return RefptsVariable(ik, "KickAngle", index=pl, name=name)

        pl = plane_(plane, "index")
        plcode = plane_(plane, "code")
        # Observables
        nm = f"{plcode}_positions"
        bpms = TrajectoryObservable(
            bpmrefs, axis=2 * pl, name=nm, target=bpmtarget, weight=bpmweight
        )
        observables = ObservableList([bpms])
        # Variables
        steeridx = ring.get_uint32_index(steerrefs)
        variables = VariableList(steerer(idx) for idx in steeridx)

        self.steerdelta = steerdelta

        super().__init__(ring, variables, observables)

    @property
    def bpmweight(self):
        """Weight of position readings."""
        return self.observables[0].weight

    @bpmweight.setter
    def bpmweight(self, value):
        self.observables[0].weight = value

    @property
    def steerdelta(self):
        """Step and weight on steerers."""
        return self.variables.deltas

    @steerdelta.setter
    def steerdelta(self, value):
        self.variables.deltas = value
