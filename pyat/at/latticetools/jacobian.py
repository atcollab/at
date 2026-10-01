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
        func: Merit function object. ``func(x)``.
        bounds: ``(n, 2)`` array of lower and upper bounds of the variables.
        max_step: Maximum change of each variable in a single step. The whole
          step is scaled down if needed, keeping its direction.
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
        max_rel_penalty_increase: Once ``n_bisections`` step-halvings
          have been tried, keep bisecting only while the penalty found
          so far is within this factor of the previous value; ``None``
          disables the extra allowance.
        broyden: If True, rank-1 (Broyden) update the Jacobian across
          steps instead of recomputing it by finite differences.
        verbose: If True, print per-step progress.
    """

    def __init__(
        self,
        func,
        bounds: np.ndarray,
        *,
        max_step: float | np.ndarray | None = None,
        n_steps_max: int = 20,
        ftol: float | None = 1e-8,
        xtol: float | None = 1e-8,
        gtol: float | None = 1e-8,
        tol: float | None = None,
        max_nfev: int | None = None,
        n_bisections: int = 3,
        max_rel_penalty_increase: float | None = 10.0,
        broyden: bool = False,
        verbose: bool = False,
    ):
        self.func = func
        self.lower, self.upper = np.asarray(bounds, dtype=float).T
        nvars = len(self.lower)
        if max_step is None:
            self.max_step = np.full(nvars, np.inf)
        else:
            self.max_step = np.broadcast_to(np.abs(max_step), (nvars,)).astype(float)
        self.n_steps_max = n_steps_max
        self.ftol = ftol
        self.xtol = xtol
        self.gtol = gtol
        self.tol = tol
        self.max_nfev = max_nfev
        self.n_bisections = n_bisections
        self.max_rel_penalty_increase = max_rel_penalty_increase
        self.broyden = broyden
        self.verbose = verbose

        self.x = None
        self._penalty_best = np.inf
        self._xbest = None

    def _evaluate(self, x):
        y = self.func(x)
        penalty_vec = self.func.penalty(y)
        penalty = float(np.sqrt(np.dot(penalty_vec, penalty_vec)))
        if penalty < self._penalty_best:
            self._penalty_best = penalty
            self._xbest = x.copy()
        return y, penalty

    def _trial_step(self, xstep, alpha):
        trial_xstep = 2.0**-alpha * xstep
        xnew = self.x - trial_xstep
        bound = np.where(xnew < self.lower, self.lower, self.upper)
        cross = (xnew < self.lower) | (xnew > self.upper)
        ratio = np.ones_like(trial_xstep)
        ratio[cross] = (self.x[cross] - bound[cross]) / trial_xstep[cross]
        scale = ratio.min(initial=1.0)
        xnew = self.x - scale * trial_xstep
        hit = cross & (ratio == scale)
        xnew[hit] = bound[hit]
        return self.x - np.clip(xnew, self.lower, self.upper)

    def solve(self, x0):
        func = self.func
        self.x = np.array(x0, dtype=float)
        if len(self.x) == 0:
            msg = "At least one variable should be present"
            raise AtError(msg)
        self._xbest = self.x.copy()
        status = 0
        jac = x_jac = y_jac = None

        y, penalty = self._evaluate(self.x)

        for istep in range(self.n_steps_max):
            if self.tol is not None and penalty < self.tol:
                status = 5
                break

            dx = None if jac is None else self.x - x_jac
            if self.broyden and dx is not None and np.dot(dx, dx) > 0:
                jac = jac + np.outer(y - y_jac - jac @ dx, dx) / np.dot(dx, dx)
            else:
                jac = func.get_jacobian(self.x)
            x_jac, y_jac = self.x.copy(), y.copy()

            grad = func.penalty(jac.T) @ func.penalty(y)
            at_lower = self.x <= self.lower
            at_upper = self.x >= self.upper
            mask = ~((at_lower & (grad > 0)) | (at_upper & (grad < 0)))
            grad[~mask] = 0.0

            if self.gtol is not None and np.linalg.norm(grad, ord=np.inf) < self.gtol:
                status = 1
                break

            while True:
                xstep = func.solve_step(jac, y, mask)
                outwards = mask & ((at_lower & (xstep > 0)) | (at_upper & (xstep < 0)))
                if not np.any(outwards):
                    break
                mask = mask & ~outwards
            with np.errstate(divide="ignore"):
                xstep = xstep * np.min(self.max_step / np.abs(xstep), initial=1.0)

            alpha = -1
            new_penalty = np.inf
            while new_penalty >= penalty:
                if alpha > self.n_bisections and (
                    self.max_rel_penalty_increase is None
                    or new_penalty < self.max_rel_penalty_increase * penalty
                ):
                    break
                alpha += 1
                trial_xstep = self._trial_step(xstep, alpha)
                y_new, new_penalty = self._evaluate(self.x - trial_xstep)

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
            y, penalty = y_new, new_penalty

            if self.verbose:
                print(f"step {istep}: alpha {alpha}, penalty {penalty}")

            if ftol_ok and xtol_ok:
                status = 4
                break
            if ftol_ok:
                status = 2
                break
            if xtol_ok:
                status = 3
                break
            if self.max_nfev is not None and func.nfev >= self.max_nfev:
                status = 0
                break

        return self._xbest, status


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
        rcond: float | None = None,
        sing_val_cutoff: int | None = None,
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
        self.rcond = rcond
        self.sing_val_cutoff = sing_val_cutoff
        self.use_mp = use_mp
        self.pool_size = pool_size
        self.start_method = start_method
        self.one_sided = one_sided
        self.nfev = 0
        self._last = None
        self._rm = ResponseMatrix(variables, constraints, ring=ring, **eval_kw)

    def close(self) -> None:
        """Shut down the multiprocessing pool, if any."""
        self._rm.close_pool()

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """Raw (unweighted) deviations from the targets at ``x``."""
        self.nfev += 1
        self.variables.set(x, ring=self._ring, **self.eval_kw)
        self.constraints.evaluate(ring=self._ring, **self.eval_kw)
        values = self.constraints.get_flat_values(err=np.nan)
        self._last = (x.copy(), values) if np.all(np.isfinite(values)) else None
        return self.constraints.get_flat_deviations(err=self.err)

    def result(self, x0: np.ndarray, x: np.ndarray) -> np.ndarray:
        """Raw deviations at ``x``, with the variables' initial values reset
        to ``x0``.
        """
        self.variables.set(x0, ring=self._ring, **self.eval_kw)
        self.variables.get(initial=True, ring=self._ring, **self.eval_kw)
        return self(x)

    def penalty(self, y: np.ndarray) -> np.ndarray:
        """Weighted deviations, used only for the solver's convergence
        checks.
        """
        return y / self.constraints.get_flat_weights()

    def solve_step(
        self, jac: np.ndarray, y: np.ndarray, mask: np.ndarray
    ) -> np.ndarray:
        """Newton step for the raw deviation ``y`` with the variables in
        ``mask``, using :meth:`.ResponseMatrix.correction_matrix`.

        ``sing_val_cutoff`` is the number of (largest) singular values to
        keep. If not given, it is derived from the relative threshold
        ``rcond``.
        """
        self._rm._varmask = mask
        self._rm.response = jac
        self._rm.solve()

        s = self._rm.singular_values
        if s.size == 0:
            return np.zeros(jac.shape[1])
        if self.sing_val_cutoff is not None:
            nvals = int(self.sing_val_cutoff)
        else:
            rcond = self.rcond
            if rcond is None:
                rcond = max(jac.shape) * np.finfo(float).eps
            nvals = int(np.sum(s > rcond * s[0]))

        return self._rm.correction_matrix(nvals=nvals) @ y

    def get_jacobian(self, x: np.ndarray) -> np.ndarray:
        """Finite-difference Jacobian at ``x``.

        :meth:`.ResponseMatrix.build` is not used because, with a pool kept
        open across Newton steps, the workers must receive the current ``x``.
        """
        self.variables.set(x, ring=self._ring, **self.eval_kw)
        self.variables.get(initial=True, ring=self._ring, **self.eval_kw)
        if self.use_mp:
            self._rm.open_pool(pool_size=self.pool_size, start_method=self.start_method)
        f0 = None
        if self.one_sided and self._last is not None:
            xlast, values = self._last
            if np.array_equal(xlast, x):
                f0 = values
        nvars = len(self.variables)
        columns = self._rm._columns(
            np.arange(nvars), x=x, one_sided=self.one_sided, f0=f0
        )
        if self.one_sided:
            self.nfev += nvars + (f0 is None)
        else:
            self.nfev += 2 * nvars
        return np.stack(columns, axis=-1)


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
        max_step: Maximum change of each variable in a single step. The whole
          step is scaled down if needed, keeping its direction.
        max_rel_penalty_increase: Once ``n_bisections`` step-halvings
          have been tried, keep bisecting only while the penalty found so
          far is within this factor of the previous value; ``None``
          disables the extra allowance.
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
       * The variable bounds are respected at every step: a step that would
         cross a bound is scaled down to stop on it, and a variable on a bound
         is held fixed as long as the gradient pushes it outwards. Near a
         bound, the finite differences are taken towards the inside.
    """
    initial_values = variables.get(initial=True, check_bounds=True, **eval_kw)

    merit = _MeritFunction(
        variables,
        constraints,
        err=err,
        rcond=rcond,
        sing_val_cutoff=sing_val_cutoff,
        use_mp=use_mp,
        pool_size=pool_size,
        start_method=start_method,
        one_sided=one_sided,
        **eval_kw,
    )
    solver = _JacobianSolver(
        merit,
        [var.bounds for var in variables],
        max_step=max_step,
        n_steps_max=n_steps_max,
        ftol=ftol,
        xtol=xtol,
        gtol=gtol,
        tol=tol,
        max_nfev=max_nfev,
        n_bisections=n_bisections,
        max_rel_penalty_increase=max_rel_penalty_increase,
        broyden=broyden,
        verbose=verbose >= 2,
    )
    try:
        xbest, status = solver.solve(initial_values)
        fbest = merit.result(initial_values, xbest)
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
