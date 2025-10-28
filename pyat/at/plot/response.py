"""Plot :py:class:`.Observable` value as a function of
:py:class:`Variable <.VariableBase>`.
"""

from __future__ import annotations

__all__ = ["plot_response"]

from collections.abc import Iterable, Mapping

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from ..lattice import VariableBase
from ..latticetools import ObservableList


def plot_response(
    var: VariableBase,
    obs: ObservableList,
    rng: Iterable[float],
    ax: Axes | None = None,
    xlabel: str = "",
    ylabel: str = "",
    title: str = "",
    color_offset: int = 0,
    **kwargs
) -> Axes:
    # noinspection PyUnresolvedReferences
    """Plot *obs* values as a function of *var*.

    Args:
        var:            Variable object.
        obs:            list of Observables.
        rng:            range of variation for the variable.
        ax:             :py:class:`~matplotlib.axes.Axes` object. If :py:obj:`None`,
          a new figure will be created.
        xlabel:         x-axis label. If empty, the variable name will be used.
        ylabel:         y-axis label.
        title:          plot title.
        color_offset:   offset in the matplotlib line color cycle.
        **kwargs:       Additional keyword arguments are transmitted to the
          :py:class:`~matplotlib.axes.Axes` creation function.

    Returns:
        ax:             the :py:class:`~matplotlib.axes.Axes` object.

    Example:
        Minimal example using only default values:

        >>> obs = at.ObservableList(
        ...     [
        ...         at.EmittanceObservable("emittances", plane="x"),
        ...         at.EmittanceObservable("emittances", plane="y"),
        ...     ],
        ...     ring=ring,
        ... )
        >>> var = at.AttributeVariable(ring, "energy", name="energy [eV]")
        >>> plot_response(var, obs, np.arange(3.0e9, 6.01e9, 0.5e9))
        >>>

        .. image:: /images/emittance_response.*
           :alt: emittance response

        Example showing the formatting possibilities by:

        - using the :py:attr:`.Observable.plot_fmt` attribute for line formatting,
        - using  the :py:attr:`.Observable.name` attribute for curve labels,
        - using dual y-axis by calling :py:func:`plot_response` twice,
        - avoiding duplicate line colors with the *color_offset* parameter,
        - using the *ylabel* and *title* parameters.

        >>> obsleft =at.ObservableList(
        ...     [
        ...         at.LocalOpticsObservable(
        ...             [0], "beta", plane="x",
        ...             name=r"$\\beta_x$",
        ...             plot_fmt={"linewidth": 3.0, "marker": "o"}
        ...         ),
        ...         at.LocalOpticsObservable(
        ...             [0], "beta", plane="y", name=r"$\\beta_z$", plot_fmt="--"
        ...         )
        ...     ],
        ...     ring=ring
        ... )
        >>>
        >>> obsright =at.ObservableList(
        ...     [
        ...         at.GlobalOpticsObservable("tune", plane="x", name=r"$\\nu_x$"),
        ...         at.GlobalOpticsObservable("tune", plane="y", name=r"$\\nu_x$"),
        ...     ],
        ...     ring=ring
        ... )
        >>> # On the left y-axis
        >>> ax = at.plot_response(
        ...     var,
        ...     obsleft,
        ...     np.arange(2.4, 2.7, 0.02),
        ...     ylabel="beta [m]",
        ...     title="Example of plot_response"
        ... )
        >>> # On the right y-axis
        >>> ax2 = at.plot_response(
        ...     var,
        ...     obsright,
        ...     np.arange(2.4, 2.7, 0.02),
        ...     ylabel="tune",
        ...     ax=ax.twinx(),
        ...     color_offset=2
        ... )
        >>> ax.set_ylim(0.0, 10.0)
        >>> ax2.set_ylim(0.0, 1.2)
        >>> ax2.grid(False)

        .. image:: /images/beta_response.*
           :alt: beta response

    """

    def compute(v):
        """Evaluate the observables for 1 variable value."""
        var.value = v
        obs.evaluate()
        return obs.values

    def plot1(x, y, obs, ncurve):
        """Plot 1 curve."""
        fmt = getattr(obs, "plot_fmt", f"C{ncurve}")
        if isinstance(fmt, Mapping):
            return ax.plot(x, y, label=obs.name, **fmt)
        else:
            return ax.plot(x, y, fmt, label=obs.name)

    if ax is None:
        _, ax = plt.subplots(subplot_kw=kwargs)
    with var.restore():
        vals = [(v, *compute(v)) for v in rng]
        xx, *yy = zip(*vals, strict=True)
        lines = [  # noqa: F841
            plot1(xx, np.array(y), ob, n + color_offset)
            for n, (y, ob) in enumerate(zip(yy, obs, strict=True))
        ]
        ax.set_xlabel(xlabel or var.name)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        ax.grid(True)
        return ax
