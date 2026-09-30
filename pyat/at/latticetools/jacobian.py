"""Jacobian (damped Gauss-Newton) matching.

The solver is adapted from the ``JacobianSolver`` of Xsuite's
``xdeps`` optimizer
(https://github.com/xsuite/xdeps/blob/main/xdeps/optimize/jacobian.py,
Apache License 2.0).

:class:`_JacobianSolver` is generic. :class:`_MeritFunction`
adapts it to AT's :class:`.VariableList` / :class:`.ObservableList`,
using :class:`.ResponseMatrix` to build the Jacobian and solve each
Newton step.
"""

from __future__ import annotations

__all__ = ["jacobian_match"]

import numpy as np
from scipy.optimize import OptimizeResult

from .observablelist import ObservableList
from .response_matrix import ResponseMatrix
from ..lattice import AtError, VariableList


_TERMINATION_MESSAGES = {
    0: "The maximum number of steps or function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`ftol` termination condition is satisfied.",
    3: "`xtol` termination condition is satisfied.",
    4: "Both `ftol` and `xtol` termination conditions are satisfied.",
    5: "`tol` termination condition is satisfied.",
    6: "No decrease of the penalty found in the line search.",
}


class _JacobianSolver:
    """Damped Gauss-Newton solver driven by a merit-function object.

    Args:
        func: Merit function object. ``func(x)`` returns the residual
          vector ``y``. It must also provide ``get_jacobian(x,
          mask_input)``, ``solve_step(...)``, ``mask_input``,
          ``mask_output``, ``_bounds`` and ``_clip_to_max_steps``. If it
          defines ``penalty(y)``, the penalty is the norm of
          ``penalty(y)`` instead of the norm of ``y``.
        n_steps_max: Maximum number of outer (Newton) steps.
        ftol: Stop when the relative decrease of the cost
          (``0.5 * penalty**2``) over a step is below this value, as in
          :func:`scipy.optimize.least_squares`. ``None`` disables it.
        xtol: Stop when ``norm(dx) < xtol * (xtol + norm(x))``, as in
          :func:`scipy.optimize.least_squares`. ``None`` disables it.
        gtol: Stop when the infinity norm of the cost gradient is below
          this value, as in :func:`scipy.optimize.least_squares`.
          ``None`` disables it.
        tol: Stop when the penalty drops below this value. ``None``
          disables it.
        max_nfev: Stop when ``func.nfev`` reaches this value. ``None``
          disables it.
        n_bisections: Maximum number of step-halvings in the line search.
        error_on_penalty_increase: Raise ``AtError`` if a line-search
          trial makes the penalty worse than this many times the
          previous value. ``False`` disables the check.
        max_rel_penalty_increase: Once ``n_bisections`` step-halvings
          have been tried, keep bisecting only while the penalty found
          so far is within this factor of the previous value; ``None``
          disables the extra allowance.
        verbose: If True, print per-step progress.
    """

    def __init__(
        self,
        func,
        n_steps_max: int = 20,
        ftol: float | None = 1e-8,
        xtol: float | None = 1e-8,
        gtol: float | None = 1e-8,
        tol: float | None = None,
        max_nfev: int | None = None,
        n_bisections: int = 3,
        error_on_penalty_increase: float | bool = 100,
        max_rel_penalty_increase: float | None = 10.0,
        verbose: bool = False,
    ):
        self.func = func
        self.n_steps_max = n_steps_max
        self.ftol = ftol
        self.xtol = xtol
        self.gtol = gtol
        self.tol = tol
        self.max_nfev = max_nfev
        self.n_bisections = n_bisections
        self.error_on_penalty_increase = error_on_penalty_increase
        self.max_rel_penalty_increase = max_rel_penalty_increase
        self.verbose = verbose

        self._x = None
        self._step = 0
        self._step_best = 0
        self._penalty_best = 1e200
        self._xbest = None
        self.mask_from_limits = None
        self._last_jac = None
        self._last_jac_x = None
        self._last_y = None

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x = np.array(np.atleast_1d(value), dtype=float)
        self.mask_from_limits = np.ones(len(self._x), dtype=bool)

    def _weighted(self, y):
        return self.func.penalty(y) if hasattr(self.func, "penalty") else y

    def eval(self, x):
        y = self.func(x)
        penalty_vec = self._weighted(y)
        penalty = float(np.sqrt(np.dot(penalty_vec, penalty_vec)))
        if self.verbose:
            print(f"penalty: {penalty}")
        if penalty < self._penalty_best:
            self._step_best = self._step
            self._penalty_best = penalty
            self._xbest = x.copy()
            if self.verbose:
                print(f"new best: {self._penalty_best}")
        return y, penalty

    def _jacobian(self, y, mask_input, broyden):
        """Broyden rank-1 update of the last Jacobian if requested and
        possible, otherwise a fresh finite-difference Jacobian.
        """
        if broyden and self._last_jac is not None:
            dx = self.x - self._last_jac_x
            dx_sqnorm = np.dot(dx, dx)
            if dx_sqnorm > 0:
                dy = y - self._last_y
                return (
                    self._last_jac + np.outer(dy - self._last_jac @ dx, dx) / dx_sqnorm
                )
        return self.func.get_jacobian(self.x, mask_input=mask_input)

    def step(
        self,
        n_steps: int = 1,
        rcond: float | None = None,
        sing_val_cutoff: int | None = None,
        broyden: bool = False,
    ):
        merit_func = self.func
        status = 0

        for step_index in range(n_steps):
            self._step += 1

            y, penalty = self.eval(self.x)

            if self.tol is not None and penalty < self.tol:
                status = 5
                break

            if len(merit_func.mask_input) == 0:
                msg = "At least one variable should be present"
                raise AtError(msg)
            if not np.any(merit_func.mask_input):
                msg = "At least one variable should be active"
                raise AtError(msg)

            mask_input = merit_func.mask_input & self.mask_from_limits

            jac = self._jacobian(y, mask_input, broyden)
            self._last_jac_x = self.x.copy()
            self._last_jac = jac.copy()
            self._last_y = y.copy()

            if self.gtol is not None:
                grad = self._weighted(jac.T) @ self._weighted(y)
                grad[~mask_input] = 0.0
                if np.linalg.norm(grad, ord=np.inf) < self.gtol:
                    status = 1
                    break

            xstep = merit_func.solve_step(
                jac,
                y,
                mask_input,
                merit_func.mask_output,
                rcond=rcond,
                sing_val_cutoff=sing_val_cutoff,
            )
            xstep = merit_func._clip_to_max_steps(xstep)
            self.mask_from_limits[:] = True

            alpha = -1
            limits = merit_func._bounds
            new_penalty = None

            while True:
                if alpha > self.n_bisections and (
                    self.max_rel_penalty_increase is None
                    or new_penalty < self.max_rel_penalty_increase * penalty
                ):
                    break
                alpha += 1
                if self.verbose:
                    print(f"\n--> step {step_index} alpha {alpha}\n")

                trial_xstep = 2.0**-alpha * xstep

                mask_hit_limit = np.zeros(len(self.x), dtype=bool)
                for ivar in range(len(self.x)):
                    xnew = self.x[ivar] - trial_xstep[ivar]
                    if xnew < limits[ivar][0]:
                        bound = limits[ivar][0]
                    elif xnew > limits[ivar][1]:
                        bound = limits[ivar][1]
                    else:
                        continue
                    trial_xstep = (
                        trial_xstep * (self.x[ivar] - bound) / trial_xstep[ivar]
                    )
                    mask_hit_limit[ivar] = True

                y, new_penalty = self.eval(self.x - trial_xstep)
                if self.verbose:
                    print(f"penalty {penalty} new_penalty {new_penalty}")

                if new_penalty < penalty:
                    break

            if (
                self.error_on_penalty_increase
                and new_penalty > penalty * self.error_on_penalty_increase
            ):
                self.eval(self.x)
                msg = (
                    f"penalty increased by more than "
                    f"{self.error_on_penalty_increase} times"
                )
                raise AtError(msg)

            if new_penalty >= penalty:
                status = 6
                break

            ftol_ok = (
                self.ftol is not None
                and penalty**2 - new_penalty**2 < self.ftol * penalty**2
            )
            xtol_ok = self.xtol is not None and np.linalg.norm(
                trial_xstep
            ) < self.xtol * (self.xtol + np.linalg.norm(self.x))

            self.x = self.x - trial_xstep
            self.mask_from_limits = ~mask_hit_limit

            if self.verbose:
                print(f"step {step_index} step_best {self._step_best} {trial_xstep}")

            if ftol_ok and xtol_ok:
                status = 4
                break
            if ftol_ok:
                status = 2
                break
            if xtol_ok:
                status = 3
                break
            nfev = getattr(merit_func, "nfev", 0)
            if self.max_nfev is not None and nfev >= self.max_nfev:
                status = 0
                break
        else:
            status = 0

        return self._xbest, status

    def solve(
        self,
        x0,
        rcond: float | None = None,
        sing_val_cutoff: int | None = None,
        broyden: bool = False,
    ):
        self.x = np.array(x0, dtype=float).copy()
        self._xbest = self.x.copy()
        return self.step(
            self.n_steps_max,
            rcond=rcond,
            sing_val_cutoff=sing_val_cutoff,
            broyden=broyden,
        )


class _MeritFunction:
    """Merit function for :class:`_JacobianSolver` built on a
    (:class:`.VariableList`, :class:`.ObservableList`) pair.

    The residual returned by ``__call__`` and the Jacobian returned by
    :meth:`get_jacobian` are both in raw (unweighted) units, as expected
    by :class:`.ResponseMatrix`, which applies the variable (``delta``)
    and observable (``weight``) weights itself when solving. The weighted
    residual used for convergence checks is given by :meth:`penalty`.
    """

    def __init__(
        self,
        variables: VariableList,
        constraints: ObservableList,
        ring,
        *,
        err: float = 1.0e6,
        max_step: np.ndarray | float | None = None,
        use_mp: bool = False,
        pool_size: int | None = None,
        start_method: str | None = None,
        one_sided: bool = True,
        **eval_kw,
    ):
        self.variables = variables
        self.constraints = constraints
        self._ring = ring
        self.eval_kw = eval_kw
        self.err = err
        self.use_mp = use_mp
        self.pool_size = pool_size
        self.start_method = start_method
        self.one_sided = one_sided
        self.nfev = 0

        n = len(variables)
        self.mask_input = np.ones(n, dtype=bool)
        self._bounds = np.array([var.bounds for var in variables], dtype=float)
        if max_step is None:
            self._max_step = np.full(n, np.inf)
        else:
            self._max_step = np.broadcast_to(
                np.abs(np.asarray(max_step, dtype=float)), (n,)
            ).copy()
        self.mask_output = None

        self._rm = ResponseMatrix(variables, constraints, ring=ring, **eval_kw)

    def close(self) -> None:
        """Shut down the multiprocessing pool, if any."""
        self._rm.close_pool()

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """Raw (unweighted) deviations from the targets at ``x``."""
        self.nfev += 1
        self.variables.set(x, ring=self._ring, **self.eval_kw)
        self.constraints.evaluate(ring=self._ring, **self.eval_kw)
        y = self.constraints.get_flat_deviations(err=self.err)
        if self.mask_output is None:
            self.mask_output = np.ones(len(y), dtype=bool)
        return y

    def penalty(self, y: np.ndarray) -> np.ndarray:
        """Weighted deviations, used only for the solver's convergence
        checks.
        """
        return y / self.constraints.get_flat_weights()

    def solve_step(
        self,
        jac: np.ndarray,
        y: np.ndarray,
        mask_input: np.ndarray,
        mask_output: np.ndarray,
        rcond: float | None = None,
        sing_val_cutoff: int | None = None,
    ) -> np.ndarray:
        """Newton step for the raw deviation ``y``, using
        :meth:`.ResponseMatrix.correction_matrix`.

        ``sing_val_cutoff`` is the number of (largest) singular values to
        keep. If not given, it is derived from the relative threshold
        ``rcond``.
        """
        nvar = jac.shape[1]
        if not np.any(mask_output):
            return np.zeros(nvar)

        self._rm._varmask = mask_input
        self._rm._obsmask = mask_output
        self._rm.response = jac
        self._rm.solve()

        s = self._rm.singular_values
        if s.size == 0:
            return np.zeros(nvar)
        if sing_val_cutoff is not None:
            nvals = int(sing_val_cutoff)
        else:
            if rcond is None:
                rcond = max(jac.shape) * np.finfo(float).eps
            nvals = int(np.sum(s > rcond * s[0]))

        return self._rm.correction_matrix(nvals=nvals) @ y

    def _clip_to_max_steps(self, xstep: np.ndarray) -> np.ndarray:
        return np.clip(xstep, -self._max_step, self._max_step)

    def get_jacobian(
        self,
        x: np.ndarray,
        mask_input: np.ndarray | None = None,
    ) -> np.ndarray:
        """Finite-difference Jacobian at ``x``.

        Only the columns of the ``mask_input``-active variables are
        computed; the others are left as zero and are excluded again at
        the SVD stage. :meth:`.ResponseMatrix.build` is not used because
        it cannot be restricted to a subset of variables.
        ``mask_input=None`` means all variables are active.
        """
        self.variables.set(x, ring=self._ring, **self.eval_kw)
        self.variables.get(initial=True, ring=self._ring, **self.eval_kw)

        nobs, nvar = self._rm.shape
        if mask_input is None:
            mask_input = np.ones(nvar, dtype=bool)
        active_idx = np.flatnonzero(mask_input)
        active_vars = [self.variables[i] for i in active_idx]

        jac = np.zeros((nobs, nvar))
        if len(active_vars) > 0:
            if self.use_mp:
                self._rm.open_pool(
                    pool_size=self.pool_size, start_method=self.start_method
                )
            columns = self._rm._columns(active_idx, x=x, one_sided=self.one_sided)
            jac[:, active_idx] = np.stack(columns, axis=-1)
            nvars = len(active_vars)
            self.nfev += nvars + 1 if self.one_sided else 2 * nvars

        return jac


def jacobian_match(
    variables: VariableList,
    constraints: ObservableList,
    *,
    verbose: int = 0,
    ftol: float | None = 1e-8,
    xtol: float | None = 1e-8,
    gtol: float | None = 1e-8,
    tol: float | None = None,
    max_nfev: int | None = None,
    n_bisections: int = 3,
    n_steps_max: int = 20,
    rcond: float | None = None,
    sing_val_cutoff: int | None = None,
    broyden: bool = False,
    max_step: float | np.ndarray | None = None,
    error_on_penalty_increase: float | bool = 100,
    max_rel_penalty_increase: float | None = 10.0,
    err: float = 1.0e6,
    use_mp: bool = False,
    pool_size: int | None = None,
    start_method: str | None = None,
    one_sided: bool = True,
    **eval_kw,
) -> OptimizeResult:
    r"""Observable matching using a damped Gauss-Newton (Jacobian) solver.

    An alternative to the scipy-backed methods of :func:`.match`, better
    suited to problems where Jacobian evaluation dominates the cost, or
    where explicit control over the singular-value cutoff or Jacobian
    reuse (``broyden=True``) is wanted.

    Args:
        variables: Variable parameters.
        constraints: Constraints to fulfill.
        verbose: 0 for silent, 1 for a summary, >=2 for per-step detail.
        ftol: Tolerance on the relative decrease of the cost
          ``0.5 * sum(weighted_deviations**2)``, as in
          :func:`scipy.optimize.least_squares`. ``None`` disables it.
        xtol: Tolerance on the relative change of the variables, as in
          :func:`scipy.optimize.least_squares`. ``None`` disables it.
        gtol: Tolerance on the infinity norm of the cost gradient, as in
          :func:`scipy.optimize.least_squares`. ``None`` disables it.
        tol: Stop when the penalty (norm of the weighted deviations)
          drops below this value. ``None`` disables it.
        max_nfev: Maximum number of function evaluations, checked after
          each Newton step. ``None`` means no limit.
        n_bisections: Maximum number of step-halvings in the line search.
        n_steps_max: Maximum number of outer (Newton) steps.
        rcond: Relative singular-value cutoff for the Newton step, used
          only when ``sing_val_cutoff`` is not given.
        sing_val_cutoff: Number of (largest) singular values to keep for
          the Newton step; takes precedence over ``rcond``.
        broyden: If True, rank-1 (Broyden) update the Jacobian across
          outer steps instead of recomputing it by finite differences.
        max_step: Absolute cap on a single Newton step, per variable.
        error_on_penalty_increase, max_rel_penalty_increase: passed
          through to :class:`_JacobianSolver`.
        err: Value substituted for a constraint's deviation when it
          cannot be evaluated.
        use_mp, pool_size, start_method: Compute the Jacobian with
          multiprocessing, as in :meth:`.ResponseMatrix.build`.
        one_sided: Compute the Jacobian using positive delta only
          (default) instead of +/-delta.

    Keyword Args:
        \*\*eval_kw: Evaluation keywords passed to
          :meth:`.VariableList.set`/`.get` and
          :meth:`.ObservableList.evaluate` (e.g. *ring*, *dp*, *dct*,
          *df*...). *ring* is required.

    Returns:
        :class:`scipy.optimize.OptimizeResult` with ``x``, ``success``,
        ``status``, ``message``, ``nfev`` and ``fun`` (the raw,
        unweighted deviations at ``x``). ``status`` follows
        :func:`scipy.optimize.least_squares` (0: maximum number of
        steps/evaluations, 1: ``gtol``, 2: ``ftol``, 3: ``xtol``, 4: both
        ``ftol`` and ``xtol``), with in addition 5: ``tol`` and 6: no
        decrease found in the line search. ``success`` is
        ``status > 0``.

    Raises:
        AtError: If the penalty increases by more than
          ``error_on_penalty_increase`` times in a single line-search
          step.

    .. note::

       * *use_mp* only parallelises the finite-difference Jacobian; the line
         search and the SVD step remain serial. Each match also pays a fixed
         cost to start the worker pool. *use_mp* pays off only when building a
         Jacobian is expensive; on small problems, it can be slower than serial.
       * With :pycode:`broyden=True` the Jacobian is computed by finite
         differences only at the first step, then updated from each accepted
         step. This saves most of the Jacobian evaluations, however on
         strongly nonlinear problems the solver may take poorer steps or
         stop early.
       * With :pycode:`broyden=True` only one Jacobian is computed, so
         :pycode:`use_mp=True` rarely helps.
    """
    initial_values = variables.get(initial=True, check_bounds=True, **eval_kw)

    merit = _MeritFunction(
        variables,
        constraints,
        err=err,
        max_step=max_step,
        use_mp=use_mp,
        pool_size=pool_size,
        start_method=start_method,
        one_sided=one_sided,
        **eval_kw,
    )
    solver = _JacobianSolver(
        merit,
        n_steps_max=n_steps_max,
        ftol=ftol,
        xtol=xtol,
        gtol=gtol,
        tol=tol,
        max_nfev=max_nfev,
        n_bisections=n_bisections,
        error_on_penalty_increase=error_on_penalty_increase,
        max_rel_penalty_increase=max_rel_penalty_increase,
        verbose=verbose >= 2,
    )
    try:
        xbest, status = solver.solve(
            initial_values,
            rcond=rcond,
            sing_val_cutoff=sing_val_cutoff,
            broyden=broyden,
        )
        variables.set(initial_values, **eval_kw)
        variables.get(initial=True, **eval_kw)
        fbest = merit(xbest)
    finally:
        merit.close()

    if verbose >= 1:
        print(f"{_TERMINATION_MESSAGES[status]} (nfev={merit.nfev})")

    return OptimizeResult(
        x=xbest,
        success=status > 0,
        status=status,
        message=_TERMINATION_MESSAGES[status],
        nfev=merit.nfev,
        fun=fbest,
    )
