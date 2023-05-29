"""Matching of lattice parameters"""
from __future__ import annotations
from ..lattice import Lattice
from ..latticetools import VariableList, ObservableList
from typing import Optional
import numpy as np
from scipy.optimize import least_squares


def match(ring: Lattice, variables: VariableList, constraints: ObservableList,
          *,
          copy: bool = True,
          method: Optional[str] = None,
          verbose: int = 2,
          max_nfev: int = 1000,
          **kwargs):
    """Perform matching of constraints by varying variables

    Minimisation is performed by the :py:func:`~scipy.optimize.least_squares`
    function.

    Args:
        ring:           Lattice description
        variables:      Variable parameters
        constraints:    Constraints to fulfill
        copy:           If :py:obj:`True`, return a modified copy of *ring*,
          otherwise modify *ring* in-place
        method:         Minimisation algorithm (see
          :py:func:`~scipy.optimize.least_squares`). If :py:obj:`None`, use
          'lm' for unbounded problems, 'trf' otherwise.
        verbose:        Level of verbosity
        max_nfev:       Maximum number of function evaluation
        **kwargs:       Keyword arguments sent to
          :py:func:`~scipy.optimize.least_squares`

    Returns:
        newring:        Modified lattice if *copy* is :py:obj:`True`

    """
    def fun(vals):
        variables.set(ring, vals)
        constraints.evaluate(ring)
        return constraints.flat_weighted_deviations

    if copy:
        boolrefs = ring.get_bool_index(None)
        for var in variables:
            boolrefs |= ring.get_bool_index(var.refpts)
        ring = ring.replace(boolrefs)

    vini = variables.get(ring, initial=True)
    bounds = np.array([var.bounds for var in variables])

    constraints.evaluate(ring, initial=True)
    ntargets = constraints.flat_values.size

    if method is None:
        if np.any(np.isfinite(bounds)) or ntargets < len(variables):
            method = 'trf'
        else:
            method = 'lm'
    if verbose >= 1:
        print('\n{} constraints, {} variables, using method {}\n'.
              format(ntargets, len(variables), method))

    least_squares(fun, vini, bounds=(bounds[:, 0], bounds[:, 1]),
                  verbose=verbose, max_nfev=max_nfev,
                  method=method, **kwargs)

    if verbose >= 1:
        print("\nConstraints:")
        print(constraints)
        print("\nVariables:")
        print(variables)

    return ring if copy else None
