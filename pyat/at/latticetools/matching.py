"""Matching of lattice parameters.

The matching acts on a :py:class:`.VariableList` to statisfy the constraints
expressed in an :py:class:`.ObservableList`.

Examples of use of such classes are available in:

* :doc:`/p/notebooks/variables`
* :doc:`/p/notebooks/observables`

Examples of matching are given in :doc:`/p/notebooks/matching`.
"""

from __future__ import annotations

__all__ = ["match"]


import numpy as np
from scipy.optimize import least_squares

from .observablelist import ObservableList
from ..lattice import Lattice, VariableList
from ..lattice.lattice_variables import ElementVariable


def _match(
    variables: VariableList,
    constraints: ObservableList,
    *,
    method: str | None = None,
    verbose: int = 2,
    max_nfev: int = 1000,
    optim_kw: dict | None = None,
    **eval_kw,
) -> np.ndarray:
    """Observable matching.

    Minimisation is performed by the :py:func:`~scipy.optimize.least_squares`
    function.

    Args:
        variables:      Variable parameters
        constraints:    Constraints to fulfill

    Keyword Args:
         method:         Minimisation algorithm (see
          :py:func:`~scipy.optimize.least_squares`). If :py:obj:`None`, use
          'lm' for unbounded problems, 'trf' otherwise,
        verbose:        Level of verbosity,
        max_nfev:       Maximum number of function evaluation,
        optim_kw:       Dictionary of optimiser keyword arguments sent to
          :py:func:`~scipy.optimize.least_squares`,
        ••eval_kw:      Evaluation keywords provided to the
          :py:meth:`ObservableList.evaluate`method. For instance *"ring"* (for
          lattice-dependent observables), *"dp"*, *"dct"*, *"orbit"*, *"twiss_in"*,
          *"r_in"*… Default values are taken from the ObservableList.

    Returns:
        Matching result: solution for the variable values.

    Raises:
        RuntimeError: If optimisation fails
    """

    def fun(vals) -> np.ndarray:
        """Evaluation function for the minimiser."""
        variables.set(vals, ring=eval_kw.get("ring"))
        constraints.evaluate(**eval_kw)
        return constraints.get_flat_weighted_deviations(err=1.0e6)

    if optim_kw is None:
        optim_kw = {}

    vini = variables.get(ring=eval_kw.get("ring"), initial=True, check_bounds=True)
    bounds = np.array([var.bounds for var in variables])
    x_scale = np.array([var.delta for var in variables])

    constraints.evaluate(initial=True, **eval_kw)
    constraints.check()
    ntargets = constraints.flat_values.size

    if method is None:
        if np.any(np.isfinite(bounds)) or ntargets < len(variables):
            method = "trf"
        else:
            method = "lm"
    if verbose >= 1:
        print(
            f"\n{ntargets} constraints, {len(variables)} variables,"
            f" using method {method}\n"
        )

    opt = least_squares(
        fun,
        vini,
        bounds=(bounds[:, 0], bounds[:, 1]),
        x_scale=x_scale,
        verbose=verbose,
        max_nfev=max_nfev,
        method=method,
        **optim_kw,
    )

    if verbose >= 1:
        print("\nConstraints:")
        print(constraints)
        print("\nVariables:")
        print(variables)

    return opt.x


def match(
    ring: Lattice,
    variables: VariableList,
    constraints: ObservableList,
    *,
    copy: bool = False,
    method: str | None = None,
    verbose: int = 2,
    max_nfev: int = 1000,
    **kwargs,
) -> Lattice | None:
    """Match constraints by varying variables.

    Minimisation is performed by the :py:func:`~scipy.optimize.least_squares`
    function.

    Args:
        ring:           Lattice description
        variables:      Variable parameters
        constraints:    Constraints to fulfill
        copy:           If :py:obj:`True`, return a modified copy of *ring*, otherwise
          perform the match in-line
        method:             Minimisation algorithm (see
          :py:func:`~scipy.optimize.least_squares`). If :py:obj:`None`, use
          'lm' for unbounded problems, 'trf' otherwise,
        verbose:            Level of verbosity,
        max_nfev:           Maximum number of function evaluation,

    Keyword Args:
        dp (float | None):  Momentum deviation. Default taken from the ObservableList,
        dct (float | None): Path lengthening. Default taken from the ObservableList,
        df (float | None):  Deviation from the nominal RF frequency. Default taken from
          the ObservableList,
        orbit (Orbit | None):  Initial orbit. Avoids looking for the closed orbit if
          it is already known. Used for :py:class:`.MatrixObservable` and
          :py:class:`.LocalOpticsObservable`. Default taken from the ObservableList,
        twiss_in:           Initial conditions for transfer line optics
          See :py:func:`.get_optics`. Used for :py:class:`.LocalOpticsObservable`.
          Default taken from the ObservableList,
        r_in (Orbit |None): Initial trajectory, used for
          :py:class:`.TrajectoryObservable`. Default taken from the ObservableList,
        **kwargs:           The other keyword arguments sent to,
          :py:func:`~scipy.optimize.least_squares`.

    Returns:
        Modified Lattice if copy=True, None otherwise

    Raises:
        TypeError: If copy=True and ElementVariable is present
        RuntimeError: If optimisation fails

    Examples:

        See :doc:`/p/notebooks/matching`
    """

    # Separate the keywords for evaluation
    eval_kw = {}
    for key in ["dp", "dct", "df", "orbit", "twiss_in", "r_in"]:
        v = kwargs.pop(key, None)
        if v is not None:
            eval_kw[key] = v

    if copy:
        # Check that there is no ElementVariable
        for var in variables:
            if isinstance(var, ElementVariable):
                msg = "When 'copy' is True, no ElementVariable is accepted."
                raise TypeError(msg)
        newring = ring.deepcopy()
        _match(
            variables,
            constraints,
            ring=newring,
            method=method,
            verbose=verbose,
            max_nfev=max_nfev,
            optim_kw=kwargs,
            **eval_kw,
        )
        return newring
    else:
        _match(
            variables,
            constraints,
            ring=ring,
            method=method,
            verbose=verbose,
            max_nfev=max_nfev,
            optim_kw=kwargs,
            **eval_kw,
        )
        return None
