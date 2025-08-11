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
    ring: Lattice,
    variables: VariableList,
    constraints: ObservableList,
    *,
    method: str | None = None,
    verbose: int = 2,
    max_nfev: int = 1000,
    dp: float | None = None,
    dct: float | None = None,
    df: float | None = None,
    **kwargs,
) -> None:
    def fun(vals):
        """Evaluation function for the minimiser."""
        variables.set(vals, ring=ring)
        constraints.evaluate(ring, dp=dp, dct=dct, df=df)
        return constraints.get_flat_weighted_deviations(err=1.0e6)

    vini = variables.get(ring=ring, initial=True, check_bounds=True)
    bounds = np.array([var.bounds for var in variables])

    constraints.evaluate(ring, initial=True, dp=dp, dct=dct, df=df)
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

    least_squares(
        fun,
        vini,
        bounds=(bounds[:, 0], bounds[:, 1]),
        verbose=verbose,
        max_nfev=max_nfev,
        method=method,
        **kwargs,
    )

    if verbose >= 1:
        print("\nConstraints:")
        print(constraints)
        print("\nVariables:")
        print(variables)


def match(
    ring: Lattice,
    variables: VariableList,
    constraints: ObservableList,
    *,
    copy: bool = False,
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

    Keyword Args:
        method:         Minimisation algorithm (see
          :py:func:`~scipy.optimize.least_squares`). If :py:obj:`None`, use
          'lm' for unbounded problems, 'trf' otherwise.
        verbose:        Level of verbosity
        max_nfev:       Maximum number of function evaluation
        dp:             Momentum deviation.
        dct:            Path lengthening.
        df:             Deviation from the nominal RF frequency.
        **kwargs:       Keyword arguments sent to
          :py:func:`~scipy.optimize.least_squares`

    Returns:
        Modified Lattice if copy=True, None otherwise

    Raises:
        TypeError: If copy=True and ElementVariable is present
        RuntimeError: If optimisation fails

    Examples:

        See :doc:`/p/notebooks/matching`
    """
    if copy:
        # Check that there is no ElementVariable
        for var in variables:
            if isinstance(var, ElementVariable):
                raise TypeError("When 'copy' is True, no ElementVariable is accepted")
        newring = ring.deepcopy()
        _match(newring, variables, constraints, **kwargs)
        return newring
    else:
        _match(ring, variables, constraints, **kwargs)
        return None
