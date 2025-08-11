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

Generic response matrix
-----------------------

The :py:class:`ResponseMatrix` class defines a general-purpose response matrix, based
on a :py:class:`.VariableList` of quantities which will be independently varied, and an
:py:class:`.ObservableList` of quantities which will be recorded for each step.

For instance let's take the horizontal displacements of all quadrupoles as variables:

>>> variables = VariableList(
...     RefptsVariable(ik, "dx", name=f"dx_{ik}", delta=0.0001)
...     for ik in ring.get_uint32_index(at.Quadrupole)
... )

The variables are the horizontal displacement ``dx`` of all quadrupoles. The variable
name is set to *dx_nnnn* where *nnnn* is the index of the quadrupole in the lattice.
The step is set to 0.0001 m.

Let's take the horizontal positions at all beam position monitors as observables:

>>> observables = at.ObservableList([at.OrbitObservable(at.Monitor, axis="x")])

This is a single observable named *orbit[x]* by default, with multiple values.

Instantiation
^^^^^^^^^^^^^

>>> resp_dx = at.ResponseMatrix(ring, variables, observables)

At that point, the response matrix is empty.

Matrix Building
^^^^^^^^^^^^^^^

The response matrix may be filled by several means:

#. Direct assignment of an array to the :py:attr:`~.ResponseMatrix.response` property.
   The shape of the array is checked.
#. :py:meth:`~ResponseMatrix.load` loads data from a file containing previously
   saved values or experimentally measured values,
#. :py:meth:`~ResponseMatrix.build_tracking` computes the matrix using tracking,
#. For some specialized response matrices a
   :py:meth:`~OrbitResponseMatrix.build_analytical` method is available.

Matrix normalisation
^^^^^^^^^^^^^^^^^^^^

To be correctly inverted, the response matrix must be correctly normalised: the norms
of its columns must be of the same order of magnitude, and similarly for the rows.

Normalisation is done by adjusting the weights :math:`w_v` for the variables
:math:`\mathbf{V}` and :math:`w_o` for the observables :math:`\mathbf{O}`.
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

Matrix pseudo-inversion
^^^^^^^^^^^^^^^^^^^^^^^

The :py:meth:`~ResponseMatrix.solve` method computes the singular values of the
weighted response matrix.

After solving, correction is available, for instance with

* :py:meth:`~ResponseMatrix.correction_matrix` which returns the correction matrix
  (pseudo-inverse of the response matrix),
* :py:meth:`~ResponseMatrix.get_correction` which returns a correction vector when
  given error values,
* :py:meth:`~ResponseMatrix.correct` which computes and optionally applies a correction
  for the provided :py:class:`.Lattice`.

Exclusion of variables and observables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Variables may be added to a set of excluded values, and similarly for observables.
Excluding an item does not change the response matrix. The values are excluded from the
pseudo-inversion of the response, possibly reducing the number of singular values.
After inversion the correction matrix is expanded to its original size by inserting
zero lines and columns at the location of excluded items. This way:

- error and correction vectors keep the same size independently of excluded values,
- excluded error values are ignored,
- excluded corrections are set to zero.

Variables can be added to the set of excluded variables using
:py:meth:`~.ResponseMatrix.exclude_vars` and observables using
:py:meth:`~.ResponseMatrix.exclude_obs`.

After excluding items, the pseudo-inverse is discarded so one must recompute it again
by calling :py:meth:`~ResponseMatrix.solve`.

The exclusion masks can be reset with  :py:meth:`~.ResponseMatrix.reset_vars` and
:py:meth:`~.ResponseMatrix.reset_obs`.
"""

from __future__ import annotations

__all__ = [
    "sequence_split",
    "ResponseMatrix",
    "OrbitResponseMatrix",
    "TrajectoryResponseMatrix",
]

import os
import multiprocessing
import concurrent.futures
import abc
import warnings
from collections.abc import Sequence, Generator, Callable
from typing import Any, ClassVar
from itertools import chain
from functools import partial
import math

import numpy as np
import numpy.typing as npt

from .observables import ElementObservable
from .observables import TrajectoryObservable, OrbitObservable, LatticeObservable
from .observables import LocalOpticsObservable, GlobalOpticsObservable
from .observablelist import ObservableList
from ..lattice import AtError, AtWarning, Refpts, Uint32Refpts, All
from ..lattice import AxisDef, plane_, Lattice, Monitor, checkattr
from ..lattice.lattice_variables import RefptsVariable
from ..lattice.variables import VariableList

FloatArray = npt.NDArray[np.float64]

_orbit_correctors = checkattr("KickAngle")

_globring: Lattice | None = None
_globobs: ObservableList | None = None

warnings.filterwarnings("always", category=AtWarning, module=__name__)


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

    _shape: tuple[int, int]
    _obsmask: npt.NDArray[bool]
    _varmask: npt.NDArray[bool]
    _response: FloatArray | None = None
    _v: FloatArray | None = None
    _uh: FloatArray | None = None
    #: Singular values of the response matrix
    singular_values: FloatArray | None = None

    def __init__(self, nobs: int, nvar: int):
        self._shape = (nobs, nvar)
        self._obsmask = np.ones(nobs, dtype=bool)
        self._varmask = np.ones(nvar, dtype=bool)

    def reset_vars(self):
        """Reset the variable exclusion mask: enable all variables"""
        self._varmask = np.ones(self.shape[1], dtype=bool)
        self._v = None
        self._uh = None
        self.singular_values = None

    def reset_obs(self):
        """Reset the observable exclusion mask: enable all observables"""
        self._obsmask = np.ones(self.shape[0], dtype=bool)
        self._v = None
        self._uh = None
        self.singular_values = None

    @property
    @abc.abstractmethod
    def varweights(self) -> np.ndarray: ...

    @property
    @abc.abstractmethod
    def obsweights(self) -> np.ndarray: ...

    @property
    def shape(self) -> tuple[int, int]:
        """Shape of the response matrix."""
        return self._shape

    def solve(self) -> None:
        """Compute the singular values of the response matrix."""
        resp = self.weighted_response
        selected = np.ix_(self._obsmask, self._varmask)
        u, s, vh = np.linalg.svd(resp[selected], full_matrices=False)
        self._v = vh.T * (1.0 / s) * self.varweights[self._varmask].reshape(-1, 1)
        self._uh = u.T / self.obsweights[self._obsmask]
        self.singular_values = s

    def check_norm(self) -> tuple[FloatArray, FloatArray]:
        """Display the norm of the rows and columns of the weighted response matrix.

        Adjusting the variables and observable weights to equalize the norms
        of rows and columns is important.

        Returns:
            obs_norms:      Norms of observables (rows)
            var_norms:      Norms of Variables (columns)
        """
        resp = self.weighted_response
        obs = np.linalg.norm(resp, axis=1)
        var = np.linalg.norm(resp, axis=0)
        print(f"max/min Observables: {np.amax(obs) / np.amin(obs)}")
        print(f"max/min Variables: {np.amax(var) / np.amin(var)}")
        return obs, var

    @property
    def response(self) -> FloatArray:
        """Response matrix."""
        resp = self._response
        if resp is None:
            raise AtError("No matrix yet: run build() or load() first")
        return resp

    @response.setter
    def response(self, response: FloatArray) -> None:
        l1, c1 = self._shape
        l2, c2 = response.shape
        if l1 != l1 or c1 != c2:
            raise ValueError(
                f"Input matrix has incompatible shape. Expected: {self.shape}"
            )
        self._response = response

    @property
    def weighted_response(self) -> FloatArray:
        """Weighted response matrix."""
        return self.response * (self.varweights / self.obsweights.reshape(-1, 1))

    def correction_matrix(self, nvals: int | None = None) -> FloatArray:
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
        cormat = np.zeros(self._shape[::-1])
        selected = np.ix_(self._varmask, self._obsmask)
        cormat[selected] = self._v[:, :nvals] @ self._uh[:nvals, :]
        return cormat

    def get_correction(
        self, observed: FloatArray, nvals: int | None = None
    ) -> FloatArray:
        """Compute the correction of the given observation.

        Args:
            observed:   Vector of observed deviations,
            nvals:      Desired number of singular values. If :py:obj:`None`, use
              all singular values

        Returns:
            corr:       Correction vector
        """
        return -self.correction_matrix(nvals=nvals) @ observed

    def save(self, file) -> None:
        """Save a response matrix in the NumPy .npy format.

        Args:
            file:   file-like object, string, or :py:class:`pathlib.Path`: File to
              which the data is saved. If file is a file-object, it must be opened in
              binary mode. If file is a string or Path, a .npy extension will
              be appended to the filename if it does not already have one.
        """
        if self._response is None:
            raise AtError("No response matrix: run build_tracking() or load() first")
        np.save(file, self._response)

    def load(self, file) -> None:
        """Load a response matrix saved in the NumPy .npy format.

        Args:
            file:   file-like object, string, or :py:class:`pathlib.Path`: the file to
              read. A file object must always be opened in binary mode.
        """
        self.response = np.load(file)


class ResponseMatrix(_SvdSolver):
    r"""Base class for response matrices.

    It is defined by any arbitrary set of :py:class:`~.variables.VariableBase`\ s and
    :py:class:`.Observable`\s

    Addition is defined on :py:class:`ResponseMatrix` objects as the addition
    of their :py:class:`~.variables.VariableBase`\ s and :py:class:`.Observable`\s to
    produce combined responses.
    """

    ring: Lattice
    variables: VariableList  #: List of matrix :py:class:`Variable <.VariableBase>`\ s
    observables: ObservableList  #: List of matrix :py:class:`.Observable`\s
    _eval_args: dict[str, Any] = {}

    def __init__(
        self,
        ring: Lattice,
        variables: VariableList,
        observables: ObservableList,
    ):
        r"""
        Args:
            ring:           Design lattice, used to compute the response
            variables:      List of :py:class:`Variable <.VariableBase>`\ s
            observables:    List of :py:class:`.Observable`\s
        """

        def limits(obslist):
            beg = 0
            for obs in obslist:
                end = beg + obs.value.size
                yield beg, end
                beg = end

        # for efficiency of parallel computation, the variable's refpts must be integer
        for var in variables:
            var.refpts = ring.get_uint32_index(var.refpts)
        self.ring = ring
        self.variables = variables
        self.observables = observables
        variables.get(ring=ring, initial=True)
        observables.evaluate(ring=ring, initial=True)
        super().__init__(len(observables.flat_values), len(variables))
        self._ob = [self._obsmask[beg:end] for beg, end in limits(self.observables)]

    def __add__(self, other: ResponseMatrix):
        if not isinstance(other, ResponseMatrix):
            raise TypeError(
                f"Cannot add {type(other).__name__} and {type(self).__name__}"
            )
        return ResponseMatrix(
            self.ring,
            VariableList(self.variables + other.variables),
            self.observables + other.observables,
        )

    def __str__(self):
        no, nv = self.shape
        return f"{type(self).__name__}({no} observables, {nv} variables)"

    @property
    def varweights(self) -> np.ndarray:
        """Variable weights."""
        return self.variables.deltas

    @property
    def obsweights(self) -> np.ndarray:
        """Observable weights."""
        return self.observables.flat_weights

    def correct(
        self, ring: Lattice, nvals: int = None, niter: int = 1, apply: bool = False
    ) -> FloatArray:
        """Compute and optionally apply the correction.

        Args:
            ring:       Lattice description. The response matrix observables
              will be evaluated for *ring* and the deviation from target will
              be corrected
            apply:      If :py:obj:`True`, apply the correction to *ring*
            niter:      Number of iterations. For more than one iteration,
              *apply* must be :py:obj:`True`
            nvals:      Desired number of singular values. If :py:obj:`None`,
              use all singular values. *nvals* may be a scalar or an iterable with
              *niter* values.

        Returns:
            correction: Vector of correction values
        """
        if niter > 1 and not apply:
            raise ValueError("If niter > 1, 'apply' must be True")
        obs = self.observables
        if apply:
            self.variables.get(ring=ring, initial=True)
        sumcorr = np.array([0.0])
        for it, nv in zip(range(niter), np.broadcast_to(nvals, (niter,))):
            print(f'step {it+1}, nvals = {nv}')
            obs.evaluate(ring, **self._eval_args)
            err = obs.flat_deviations
            if np.any(np.isnan(err)):
                raise AtError(
                    f"Step {it + 1}: Invalid observables, cannot compute correction"
                )
            corr = self.get_correction(obs.flat_deviations, nvals=nv)
            sumcorr = sumcorr + corr  # non-broadcastable sumcorr
            if apply:
                self.variables.increment(corr, ring=ring)
        return sumcorr

    def build_tracking(
        self,
        use_mp: bool = False,
        pool_size: int | None = None,
        start_method: str | None = None,
        **kwargs,
    ) -> FloatArray:
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

        Returns:
            response:       Response matrix
        """
        self._eval_args = kwargs
        self.observables.evaluate(self.ring)
        ring = self.ring.deepcopy()

        if use_mp:
            global _globring
            global _globobs
            ctx = multiprocessing.get_context(start_method)
            if pool_size is None:
                pool_size = min(len(self.variables), os.cpu_count())
            obschunks = sequence_split(self.variables, pool_size)
            if ctx.get_start_method() == "fork":
                _globring = ring
                _globobs = self.observables
                _single_resp = partial(_resp_fork, **kwargs)
            else:
                _single_resp = partial(_resp, ring, self.observables, **kwargs)
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=pool_size,
                mp_context=ctx,
            ) as pool:
                results = list(chain(*pool.map(_single_resp, obschunks)))
            _globring = None
            _globobs = None
        else:
            results = _resp(ring, self.observables, self.variables, **kwargs)

        resp = np.stack(results, axis=-1)
        self.response = resp
        return resp

    def build_analytical(self) -> FloatArray:
        """Build the response matrix."""
        raise NotImplementedError(
            f"build_analytical not implemented for {self.__class__.__name__}"
        )

    def _on_obs(self, fun: Callable, *args, obsid: int | str = 0):
        """Apply a function to the selected observable"""
        if not isinstance(obsid, str):
            return fun(self.observables[obsid], *args)
        else:
            for obs in self.observables:
                if obs.name == obsid:
                    return fun(obs, *args)
            else:
                raise ValueError(f"Observable {obsid} not found")

    def get_target(self, *, obsid: int | str = 0) -> FloatArray:
        r"""Return the target of the specified observable

        Args:
            obsid:  :py:class:`.Observable` name or index in the observable list.

        Returns:
            target: observable target
        """
        def _get(obs):
            return obs.target

        return self._on_obs(_get, obsid=obsid)

    def set_target(self, target: npt.ArrayLike, *, obsid: int | str = 0) -> None:
        r"""Set the target of the specified observable

        Args:
            target: observable target. Must be broadcastable to the shape of the
              observable value.
            obsid:  :py:class:`.Observable` name or index in the observable list.
        """

        def _set(obs, targ):
            obs.target = targ

        return self._on_obs(_set, target, obsid=obsid)

    def exclude_obs(self, *, obsid: int | str = 0, refpts: Refpts = None) -> None:
        # noinspection PyUnresolvedReferences
        r"""Add an observable item to the set of excluded values

        After excluding observation points, the matrix must be inverted again using
        :py:meth:`solve`.

        Args:
            obsid:      :py:class:`.Observable` name or index in the observable list.
            refpts:     location of elements to exclude for
              :py:class:`.ElementObservable` objects, otherwise ignored.

        Raises:
            ValueError: No observable with the given name.
            IndexError: Observableindex out of range.

        Example:
            >>> resp = OrbitResponseMatrix(ring, "h", Monitor, Corrector)
            >>> resp.exclude_obs(obsid="x_orbit", refpts="BPM_02")

            Create an horizontal :py:class:`OrbitResponseMatrix` from
            :py:class:`.Corrector` elements to :py:class:`.Monitor` elements,
            and exclude the monitor with name "BPM_02"
        """

        def exclude(ob, msk):
            inimask = msk.copy()
            if isinstance(ob, ElementObservable) and not ob.summary:
                boolref = self.ring.get_bool_index(refpts)
                # noinspection PyProtectedMember
                msk &= np.logical_not(boolref[ob._boolrefs])
            else:
                msk[:] = False
            if np.all(msk == inimask):
                warnings.warn(AtWarning("No new excluded value"), stacklevel=3)
            # Force a new computation
            self.singular_values = None

        if not isinstance(obsid, str):
            exclude(self.observables[obsid], self._ob[obsid])
        else:
            for obs, mask in zip(self.observables, self._ob):
                if obs.name == obsid:
                    exclude(obs, mask)
                    break
            else:
                raise ValueError(f"Observable {obsid} not found")

    @property
    def excluded_obs(self) -> dict:
        """Directory of excluded observables.

        The dictionary keys are the observable names, the values are the integer
        indices of excluded items (empty list if no exclusion).
        """

        def ex(obs, mask):
            if isinstance(obs, ElementObservable) and not obs.summary:
                refpts = self.ring.get_bool_index(None)
                # noinspection PyProtectedMember
                refpts[obs._boolrefs] = np.logical_not(mask)
                refpts = self.ring.get_uint32_index(refpts)
            else:
                refpts = np.arange(0 if np.all(mask) else mask.size, dtype=np.uint32)
            return refpts

        return {ob.name: ex(ob, mask) for ob, mask in zip(self.observables, self._ob)}

    def exclude_vars(self, *varid: int | str) -> None:
        # noinspection PyUnresolvedReferences
        """Add variables to the set of excluded variables.

        Args:
            *varid:  :py:class:`Variable <.VariableBase>` names or variable indices
              in the variable list

        After excluding variables, the matrix must be inverted again using
        :py:meth:`solve`.

        Examples:
            >>> resp.exclude_vars(0, "var1", -1)

            Exclude the 1st variable, the variable named "var1" and the last variable.
        """
        nameset = set(nm for nm in varid if isinstance(nm, str))
        varidx = [nm for nm in varid if not isinstance(nm, str)]
        mask = np.array([var.name in nameset for var in self.variables])
        mask[varidx] = True
        miss = nameset - {var.name for var, ok in zip(self.variables, mask) if ok}
        if miss:
            raise ValueError(f"Unknown variables: {miss}")
        self._varmask &= np.logical_not(mask)

    @property
    def excluded_vars(self) -> list:
        """List of excluded variables"""
        return [var.name for var, ok in zip(self.variables, self._varmask) if not ok]


class OrbitResponseMatrix(ResponseMatrix):
    # noinspection PyUnresolvedReferences
    r"""Orbit response matrix.

    An :py:class:`OrbitResponseMatrix` applies to a single plane, horizontal or
    vertical. A combined response matrix is obtained by adding horizontal and
    vertical matrices. However, the resulting matrix has the :py:class:`ResponseMatrix`
    class, which implies that the :py:class:`OrbitResponseMatrix` specific methods are
    not available.

    Variables are a set of steerers and optionally the RF frequency. Steerer
    variables are named ``xnnnn`` or ``ynnnn`` where nnnn is the index in the
    lattice. The RF frequency variable is named ``RF frequency``.

    Observables are the closed orbit position at selected points, named
    ``orbit[x]`` for the horizontal plane or ``orbit[y]`` for the vertical plane,
    and optionally the sum of steerer angles named ``sum(h_kicks)`` or
    ``sum(v_kicks)``

    The variable elements must have the *KickAngle* attribute used for correction.
    It's available for all magnets, though not present by default
    except in :py:class:`.Corrector` magnets. For other magnets, the attribute
    should be explicitly created.

    By default, the observables are all the :py:class:`.Monitor` elements, and the
    variables are all the elements having a *KickAngle* attribute.
    This is equivalent to:

    >>> resp_v = OrbitResponseMatrix(
    ...     ring, "v", bpmrefs=at.Monitor, steerrefs=at.checkattr("KickAngle")
    ... )
    """

    bpmrefs: Uint32Refpts  #: location of position monitors
    steerrefs: Uint32Refpts  #: location of steerers

    def __init__(
        self,
        ring: Lattice,
        plane: AxisDef,
        bpmrefs: Refpts = Monitor,
        steerrefs: Refpts = _orbit_correctors,
        *,
        cavrefs: Refpts = None,
        bpmweight: float | Sequence[float] = 1.0,
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
              is also the cavity weight. Default: automatically computed.
            steerdelta: Step on steerers for matrix computation [rad]. This is
              also the steerer weight. Must be broadcastable to the number of steerers.
            steersum:   If :py:obj:`True`, the sum of steerers is appended to the
              Observables.
            stsumweight: Weight on steerer summation. Default: automatically computed.

        :ivar VariableList variables: matrix variables
        :ivar ObservableList observables: matrix observables

        By default, the weights of cavities and steerers summation are set to give
        a factor 2 more efficiency than steerers and BPMs

        """

        def steerer(ik, delta):
            name = f"{plcode}{ik:04}"
            return RefptsVariable(ik, "KickAngle", index=pl, name=name, delta=delta)

        def set_norm():
            bpm = LocalOpticsObservable(bpmrefs, "beta", plane=pl)
            sts = LocalOpticsObservable(steerrefs, "beta", plane=pl)
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
            cd = vv * korb * alpha * freq / np.linalg.norm(result[2] / bpmweight)
            sw = np.linalg.norm(deltas) / vo / korb
            return cd, sw

        pl = plane_(plane, key="index")
        plcode = plane_(plane, key="code")
        ids = ring.get_uint32_index(steerrefs)
        nbsteers = len(ids)
        deltas = np.broadcast_to(steerdelta, nbsteers)
        if steersum and stsumweight is None or cavrefs and cavdelta is None:
            cavd, stsw = set_norm()

        # Observables
        bpms = OrbitObservable(bpmrefs, axis=2 * pl, target=bpmtarget, weight=bpmweight)
        observables = ObservableList([bpms])
        if steersum:
            # noinspection PyUnboundLocalVariable
            sumobs = LatticeObservable(
                steerrefs,
                "KickAngle",
                name=f"{plcode}_kicks",
                target=0.0,
                index=pl,
                weight=stsumweight if stsumweight else stsw / 2.0,
                statfun=np.sum,
            )
            observables.append(sumobs)

        # Variables
        variables = VariableList(steerer(ik, delta) for ik, delta in zip(ids, deltas))
        if cavrefs is not None:
            active = (el.longt_motion for el in ring.select(cavrefs))
            if not all(active):
                raise ValueError("Cavities are not active")
            # noinspection PyUnboundLocalVariable
            cavvar = RefptsVariable(
                cavrefs,
                "Frequency",
                name="RF frequency",
                delta=cavdelta if cavdelta else 2.0 * cavd,
            )
            variables.append(cavvar)

        super().__init__(ring, variables, observables)
        self.plane = pl
        self.steerrefs = ids
        self.nbsteers = nbsteers
        self.bpmrefs = ring.get_uint32_index(bpmrefs)

    def exclude_obs(self, *, obsid: int | str = 0, refpts: Refpts = None) -> None:
        # noinspection PyUnresolvedReferences
        r"""Add an observable item to the set of excluded values.

        After excluding observation points, the matrix must be inverted again using
        :py:meth:`solve`.

        Args:
            obsid:      If 0 (default), act on Monitors. Otherwise,
              it must be 1 or "sum(x_kicks)"  or "sum(y_kicks)"
            refpts:    location of Monitors to exclude

        Raises:
            ValueError: No observable with the given name.
            IndexError: Observableindex out of range.

        Example:
            >>> resp = OrbitResponseMatrix(ring, "h")
            >>> resp.exclude_obs("BPM_02")

            Create an horizontal :py:class:`OrbitResponseMatrix` from
            :py:class:`.Corrector` elements to :py:class:`.Monitor` elements,
            and exclude all monitors with name "BPM_02"
        """
        super().exclude_obs(obsid=obsid, refpts=refpts)

    def exclude_vars(self, *varid: int | str, refpts: Refpts = None) -> None:
        # noinspection PyUnresolvedReferences
        """Add correctors to the set of excluded variables.

        Args:
            *varid:  :py:class:`Variable <.VariableBase>` names or variable indices
              in the variable list
            refpts:  location of correctors to exclude

        After excluding correctors, the matrix must be inverted again using
        :py:meth:`solve`.

        Examples:
            >>> resp.exclude_vars(0, "x0097", -1)

            Exclude the 1st variable, the variable named "x0097" and the last variable.

            >>> resp.exclude_vars(refpts="SD1E")

            Exclude all variables associated with the element named "SD1E".
        """
        plcode = plane_(self.plane, key="code")
        names = [f"{plcode}{ik:04}" for ik in self.ring.get_uint32_index(refpts)]
        super().exclude_vars(*varid, *names)

    def normalise(
        self, cav_ampl: float | None = 2.0, stsum_ampl: float | None = 2.0
    ) -> None:
        """Normalise the response matrix

        Adjust the RF cavity delta and/or the weight of steerer summation so that the
        weighted response matrix is normalised.

        Args:
            cav_ampl: Desired ratio between the cavity response and the average of
              steerer responses. If :py:obj:`None`, do not normalise.
            stsum_ampl: Desired inverse ratio between the weight of the steerer
              summation and the average of Monitor responses. If :py:obj:`None`,
              do not normalise.

        By default, the normalisation gives to the RF frequency and steerer summation
        a factor 2 more efficiency than steerers and BPMs
        """
        resp = self.weighted_response
        normvar = np.linalg.norm(resp, axis=0)
        normobs = np.linalg.norm(resp, axis=1)
        if len(self.variables) > self.nbsteers and cav_ampl is not None:
            self.cavdelta *= np.mean(normvar[:-1]) / normvar[-1] * cav_ampl
        if len(self.observables) > 1 and stsum_ampl is not None:
            self.stsumweight = (
                self.stsumweight * normobs[-1] / np.mean(normobs[:-1]) / stsum_ampl
            )

    def build_analytical(self, **kwargs) -> FloatArray:
        """Build analytically the response matrix.

        Keyword Args:
            dp (float):     Momentum deviation. Defaults to :py:obj:`None`
            dct (float):    Path lengthening. Defaults to :py:obj:`None`
            df (float):     Deviation from the nominal RF frequency.
              Defaults to :py:obj:`None`

        Returns:
            response:       Response matrix

        References:
            .. [#Franchi] A. Franchi, S.M. Liuzzo, Z. Marti, *"Analytic formulas for
               the rapid evaluation of the orbit response matrix and chromatic functions
               from lattice parameters in circular accelerators"*,
               arXiv:1711.06589 [physics.acc-ph]
        """

        def tauwj(muj, muw):
            tau = muj - muw
            if tau < 0.0:
                tau += 2.0 * pi_tune
            return tau - pi_tune

        ring = self.ring
        pl = self.plane
        _, ringdata, elemdata = ring.linopt6(All, **kwargs)
        pi_tune = math.pi * ringdata.tune[pl]
        dataw = elemdata[self.steerrefs]
        dataj = elemdata[self.bpmrefs]
        dispj = dataj.dispersion[:, 2 * pl]
        dispw = dataw.dispersion[:, 2 * pl]
        lw = np.array([elem.Length for elem in ring.select(self.steerrefs)])
        taufunc = np.frompyfunc(tauwj, 2, 1)

        sqbetaw = np.sqrt(dataw.beta[:, pl])
        ts = lw / sqbetaw / 2.0
        tc = sqbetaw - dataw.alpha[:, pl] * ts
        twj = np.astype(taufunc.outer(dataj.mu[:, pl], dataw.mu[:, pl]), np.float64)
        jcwj = tc * np.cos(twj) + ts * np.sin(twj)
        coefj = np.sqrt(dataj.beta[:, pl]) / (2.0 * np.sin(pi_tune))
        resp = coefj[:, np.newaxis] * jcwj
        if ring.is_6d:
            alpha_c = ring.disable_6d(copy=True).get_mcf()
            resp += np.outer(dispj, dispw) / (alpha_c * ring.circumference)
            if len(self.variables) > self.nbsteers:
                rfrsp = -dispj / (alpha_c * ring.rf_frequency)
                resp = np.concatenate((resp, rfrsp[:, np.newaxis]), axis=1)
        if len(self.observables) > 1:
            sumst = np.ones(resp.shape[1], np.float64)
            if len(self.variables) > self.nbsteers:
                sumst[-1] = 0.0
            resp = np.concatenate((resp, sumst[np.newaxis]), axis=0)
        self.response = resp
        return resp

    @property
    def bpmweight(self) -> FloatArray:
        """Weight of position readings."""
        return self.observables[0].weight

    @bpmweight.setter
    def bpmweight(self, value: npt.ArrayLike):
        self.observables[0].weight = value

    @property
    def stsumweight(self) -> FloatArray:
        """Weight of steerer summation."""
        return self.observables[1].weight

    @stsumweight.setter
    def stsumweight(self, value: float):
        self.observables[1].weight = value

    @property
    def steerdelta(self) -> FloatArray:
        """Step and weight of steerers."""
        return self.variables[: self.nbsteers].deltas

    @steerdelta.setter
    def steerdelta(self, value: npt.ArrayLike):
        self.variables[: self.nbsteers].deltas = value

    @property
    def cavdelta(self) -> FloatArray:
        """Step and weight of RF frequency deviation."""
        return self.variables[self.nbsteers].delta

    @cavdelta.setter
    def cavdelta(self, value: float):
        self.variables[self.nbsteers].delta = value


class TrajectoryResponseMatrix(ResponseMatrix):
    """Trajectory response matrix.

    A :py:class:`TrajectoryResponseMatrix` applies to a single plane, horizontal or
    vertical. A combined response matrix is obtained by adding horizontal and vertical
    matrices. However, the resulting matrix has the :py:class:`ResponseMatrix`
    class, which implies that the :py:class:`OrbitResponseMatrix` specific methods are
    not available.

    Variables are a set of steerers. Steerer variables are named ``xnnnn`` or
    ``ynnnn`` where *nnnn* is the index in the lattice.

    Observables are the trajectory position at selected points, named ``trajectory[x]``
    for the horizontal plane or ``trajectory[y]`` for the vertical plane.

    The variable elements must have the *KickAngle* attribute used for correction.
    It's available for all magnets, though not present by default
    except in :py:class:`.Corrector` magnets. For other magnets, the attribute
    should be explicitly created.

    By default, the observables are all the :py:class:`.Monitor` elements, and the
    variables are all the elements having a *KickAngle* attribute.

    """

    bpmrefs: Uint32Refpts
    steerrefs: Uint32Refpts
    _default_twiss_in: ClassVar[dict] = {"beta": np.ones(2), "alpha": np.zeros(2)}

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

        def steerer(ik, delta):
            name = f"{plcode}{ik:04}"
            return RefptsVariable(ik, "KickAngle", index=pl, name=name, delta=delta)

        pl = plane_(plane, key="index")
        plcode = plane_(plane, key="code")
        ids = ring.get_uint32_index(steerrefs)
        nbsteers = len(ids)
        deltas = np.broadcast_to(steerdelta, nbsteers)
        # Observables
        bpms = TrajectoryObservable(
            bpmrefs, axis=2 * pl, target=bpmtarget, weight=bpmweight
        )
        observables = ObservableList([bpms])
        # Variables
        variables = VariableList(steerer(ik, delta) for ik, delta in zip(ids, deltas))

        super().__init__(ring, variables, observables)
        self.plane = pl
        self.steerrefs = ids
        self.nbsteers = nbsteers
        self.bpmrefs = ring.get_uint32_index(bpmrefs)

    def build_analytical(self, **kwargs) -> FloatArray:
        """Build analytically the response matrix.

        Keyword Args:
            dp (float):     Momentum deviation. Defaults to :py:obj:`None`
            dct (float):    Path lengthening. Defaults to :py:obj:`None`
            df (float):     Deviation from the nominal RF frequency.
              Defaults to :py:obj:`None`

        Returns:
            response:       Response matrix
        """
        ring = self.ring
        pl = self.plane
        twiss_in = self._eval_args.get("twiss_in", self._default_twiss_in)
        _, _, elemdata = ring.linopt6(All, twiss_in=twiss_in, **kwargs)
        dataj = elemdata[self.bpmrefs]
        dataw = elemdata[self.steerrefs]
        lw = np.array([elem.Length for elem in ring.select(self.steerrefs)])

        sqbetaw = np.sqrt(dataw.beta[:, pl])
        ts = lw / sqbetaw / 2.0
        tc = sqbetaw - dataw.alpha[:, pl] * ts
        twj = dataj.mu[:, pl].reshape(-1, 1) - dataw.mu[:, pl]
        jswj = tc * np.sin(twj) - ts * np.cos(twj)
        coefj = np.sqrt(dataj.beta[:, pl])
        resp = coefj[:, np.newaxis] * jswj
        resp[twj < 0.0] = 0.0
        self.response = resp
        return resp

    def exclude_obs(self, *, obsid: int | str = 0, refpts: Refpts = None) -> None:
        # noinspection PyUnresolvedReferences
        r"""Add a monitor to the set of excluded values.

        After excluding observation points, the matrix must be inverted again using
        :py:meth:`solve`.

        Args:
            refpts:    location of Monitors to exclude

        Raises:
            ValueError: No observable with the given name.
            IndexError: Observableindex out of range.

        Example:
            >>> resp = TrajectoryResponseMatrix(ring, "v")
            >>> resp.exclude_obs("BPM_02")

            Create a vertical :py:class:`TrajectoryResponseMatrix` from
            :py:class:`.Corrector` elements to :py:class:`.Monitor` elements,
            and exclude all monitors with name "BPM_02"
        """
        super().exclude_obs(obsid=0, refpts=refpts)

    def exclude_vars(self, *varid: int | str, refpts: Refpts = None) -> None:
        # noinspection PyUnresolvedReferences
        """Add correctors to the set of excluded variables.

        Args:
            *varid:  :py:class:`Variable <.VariableBase>` names or variable indices
              in the variable list
            refpts:  location of correctors to exclude

        After excluding correctors, the matrix must be inverted again using
        :py:meth:`solve`.

        Examples:
            >>> resp.exclude_vars(0, "x0103", -1)

            Exclude the 1st variable, the variable named "x0103" and the last variable.

            >>> resp.exclude_vars(refpts="SD1E")

            Exclude all variables associated with the element named "SD1E".
        """
        plcode = plane_(self.plane, key="code")
        names = [f"{plcode}{ik:04}" for ik in self.ring.get_uint32_index(refpts)]
        super().exclude_vars(*varid, *names)

    @property
    def bpmweight(self) -> FloatArray:
        """Weight of position readings."""
        return self.observables[0].weight

    @bpmweight.setter
    def bpmweight(self, value: npt.ArrayLike):
        self.observables[0].weight = value

    @property
    def steerdelta(self) -> np.ndarray:
        """Step and weight on steerers."""
        return self.variables.deltas

    @steerdelta.setter
    def steerdelta(self, value):
        self.variables.deltas = value
