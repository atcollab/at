"""Matching of lattice parameters.

The matching acts on a :py:class:`.VariableList` to statisfy the constraints
expressed in an :py:class:`.ObservableList`.

Examples of use of such classes are available in:

* :doc:`/p/notebooks/variables`
* :doc:`/p/notebooks/observables`

Examples of matching are given in :doc:`/p/notebooks/matching`.
"""

from __future__ import annotations

__all__ = ["MatchingMethod", "match", "ring_match"]

import inspect
from enum import Enum

import numpy as np

from .optimizers import jacobian_match, least_squares_match
from .observablelist import ObservableList
from ..lattice import Lattice, VariableList, AtError
from ..lattice.lattice_variables import ElementVariable


class MatchingMethod(Enum):
    """Methods for the matching algorithm.

    * ``LEASTSQUARES``: :py:func:`.least_squares_match`, based on
      :py:func:`~scipy.optimize.least_squares`,
    * ``JACOBIAN``: :py:func:`.jacobian_match`.
    """

    LEASTSQUARES = 0
    JACOBIAN = 1


def match(
    variables: VariableList,
    constraints: ObservableList,
    *,
    method: MatchingMethod = MatchingMethod.LEASTSQUARES,
    verbose: int = 2,
    max_nfev: int = 1000,
    lstsq_method: str | None = None,
    optim_kw: dict | None = None,
    **eval_kw,
) -> np.ndarray:
    r"""Observable matching.

    Minimisation is performed by :py:func:`.least_squares_match` or
    :py:func:`.jacobian_match`, depending on *method*.

    Args:
        variables:      Variable parameters
        constraints:    Constraints to fulfill
        method:         Matching algorithm: :py:attr:`MatchingMethod.LEASTSQUARES` for
          :py:func:`.least_squares_match` or :py:attr:`MatchingMethod.JACOBIAN`
          for :py:func:`.jacobian_match`,
        verbose:        Level of verbosity,
        max_nfev:       Maximum number of function evaluation,
        lstsq_method: :py:func:`~scipy.optimize.least_squares` algorithm,
          'trf', 'dogbox' or 'lm'. If :py:obj:`None`, use 'lm' for unbounded
          problems with at least as many constraints as variables, 'trf'
          otherwise. 'lm' supports neither bounds nor fewer constraints than
          variables.
        optim_kw:       Dictionary of optimiser keyword arguments sent to
          :py:func:`.least_squares_match` or :py:func:`.jacobian_match`. The
          tolerances *ftol*, *xtol* and *gtol* have the same meaning for both
          methods,

    Keyword Args:
        \*\*eval_kw:    Evaluation keywords provided to the
          :py:meth:`ObservableList.evaluate` method. For instance *"ring"* (for
          lattice-dependent observables), *"dp"*, *"dct"*, *"orbit"*, *"twiss_in"*,
          *"r_in"*… Default values are taken from the ObservableList.

    Returns:
        Matching result: solution for the variable values.

    Raises:
        RuntimeError: If optimisation fails
    """

    if np.any(np.isnan(constraints.get_flat_targets())):
        msg = (
            "The targets are undefined for some constraints. "
            "The residual cannot be calculated. "
            "Please set a value (float or None)"
        )
        raise AtError(msg)

    if optim_kw is None:
        optim_kw = {}

    constraints.evaluate(initial=True, **eval_kw)
    constraints.check()
    ntargets = constraints.flat_values.size

    if verbose >= 1:
        msg = (
            f"\n{ntargets} constraints, {len(variables)} variables,"
            f" using method {method.name}"
        )
        print(msg + "\n")

    if method == MatchingMethod.JACOBIAN:
        solver = jacobian_match
    elif method == MatchingMethod.LEASTSQUARES:
        solver = least_squares_match
        optim_kw = {"lstsq_method": lstsq_method, **optim_kw}
    else:
        msg = f"Undefined matching method: {method}"
        raise AtError(msg)
    unknown = set(optim_kw) - set(inspect.signature(solver).parameters)
    if unknown:
        msg = (
            f"Unexpected keyword(s) for {solver.__name__}: {', '.join(sorted(unknown))}"
        )
        raise TypeError(msg)
    opt = solver(
        variables,
        constraints,
        verbose=verbose,
        max_nfev=max_nfev,
        **optim_kw,
        **eval_kw,
    )

    if verbose >= 1:
        print("\nConstraints:")
        print(constraints)
        print("\nVariables:")
        print(variables)

    return opt.x


def ring_match(
    ring: Lattice,
    variables: VariableList,
    constraints: ObservableList,
    *,
    copy: bool = False,
    method: MatchingMethod = MatchingMethod.LEASTSQUARES,
    verbose: int = 2,
    max_nfev: int = 1000,
    lstsq_method: str | None = None,
    **kwargs,
) -> Lattice | None:
    r"""Match constraints by varying variables.

    Minimisation is performed by :py:func:`.least_squares_match` or
    :py:func:`.jacobian_match`, depending on *method*.

    Args:
        ring:           Lattice description
        variables:      Variable parameters
        constraints:    Constraints to fulfill
        copy:           If :py:obj:`True`, return a modified copy of *ring*, otherwise
          perform the match in-line
        method:             Matching algorithm: :py:attr:`MatchingMethod.LEASTSQUARES`
          or :py:attr:`MatchingMethod.JACOBIAN`, see :py:func:`match`,
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
        r_in (Orbit | None): Initial trajectory, used for
          :py:class:`.TrajectoryObservable`. Default taken from the ObservableList,
        lstsq_method: :py:func:`~scipy.optimize.least_squares` algorithm,
          'trf', 'dogbox' or 'lm'. If :py:obj:`None`, use 'lm' for unbounded
          problems with at least as many constraints as variables, 'trf'
          otherwise. 'lm' supports neither bounds nor fewer constraints than
          variables.
        \*\*kwargs:         The other keyword arguments are sent to the
          matching algorithm, see *optim_kw* in :py:func:`match`.

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
        match(
            variables,
            constraints,
            ring=newring,
            method=method,
            verbose=verbose,
            max_nfev=max_nfev,
            lstsq_method=lstsq_method,
            optim_kw=kwargs,
            **eval_kw,
        )
        return newring
    else:
        match(
            variables,
            constraints,
            ring=ring,
            method=method,
            verbose=verbose,
            max_nfev=max_nfev,
            lstsq_method=lstsq_method,
            optim_kw=kwargs,
            **eval_kw,
        )
        return None


Lattice.match = ring_match
