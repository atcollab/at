"""Plot :py:class:`.Observable` value as a function of
:py:class:`Variable <.VariableBase>`.
"""

from __future__ import annotations

__all__ = ["plot_response"]

from collections.abc import Generator, Iterable
from contextlib import contextmanager

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from ..lattice import VariableBase
from ..latticetools import ObservableList


@contextmanager
def _save_variables(var: VariableBase) -> Generator[None, None, None]:
    print("Saving variable")
    var.get(initial=True)
    try:
        yield
    finally:
        print("Restoring variable")
        var.reset()


def plot_response(
    var: VariableBase,
    obs: ObservableList,
    rng: Iterable[float],
    ax: Axes | None = None,
    xlabel: str = "",
    ylabel: str = "",
    title: str = "",
):
    """Plot *obs* values as a function of *var*.

    Args:
        var: Variable object.
        obs: list of Observables.
        rng: range of variation for the variable.
        ax: :py:class:`~matplotlib.axes.Axes` object. If :py:obj:`None`, a new figure
          will be created.
        xlabel: x-axis label. If empty, the variable name will be used.
        ylabel: y-axis label.
        title: plot title.

    Example:
        >>> obs = at.ObservableList(
        ...    [at.EmittanceObservable("emittances", plane="x"),
        ...     at.EmittanceObservable("emittances", plane="y")],
        ...    ring=ring
        ... )
        >>> var = at.AttributeVariable(ring, "energy", name="energy [eV]")
        >>> plot_response(
        ...     var, obs, np.arange(3.0e9, 6.01e9, 0.5e9),
        ...     ylabel="Emittance [m]")
        ... )

        .. image:: /images/emittance_response.*
           :alt: emittance response

        >>> obs = at.ObservableList(
        ...    [at.LocalOpticsObservable([0], "beta", plane="x"),
        ...     at.LocalOpticsObservable([0], "beta", plane="y")],
        ...    ring=ring
        ... )
        >>> var = at.RefptsVariable(
        ...     "QF1[AE]", "PolynomB", index=1, name = "QF1 strength", ring=ring
        ... )
        >>> at.plot_response(var, obs, np.arange(2.4, 2.7, 0.02), ylabel="beta[m]")

        .. image:: /images/beta_response.*
           :alt: beta response

    """

    def compute(v):
        var.value = v
        obs.evaluate()
        return obs.values

    if ax is None:
        _, ax = plt.subplots()
    with var.restore():
        vals = [(v, *compute(v)) for v in rng]
        xx, *yy = zip(*vals, strict=True)
        for hy, ob in zip(yy, obs, strict=True):
            ax.plot(xx, np.array(hy), label=ob.name)
        ax.set_xlabel(xlabel or var.name)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        ax.grid(True)
