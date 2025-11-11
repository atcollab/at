"""Plot :py:class:`.Observable` values as a function of a
:py:class:`Variable <.VariableBase>`.
"""

from __future__ import annotations

__all__ = ["plot_response"]

from collections.abc import Mapping
import itertools

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from ..lattice import VariableBase
from ..latticetools import ObservableList


def plot_response(
    var: VariableBase,
    *args,
    axes: Axes | None = None,
    xlabel: str = "",
    ylabel: str = "",
    **kwargs
) -> Axes:
    # noinspection PyUnresolvedReferences
    r"""plot_response(var: VariableBase, rng: Iterable[float], obsleft: ObservableList, obsright: ObservableList, **kwargs) -> Axes
    Plot :py:class:`.Observable` values as a function of a :py:class:`Variable <.VariableBase>`.

    Args:
        var:        Variable object,
        rng:        range of variation for the variable,
        obsleft:    List of Observables plotted on the left axis,
        obsright:   Optional list of Observables plotted on the right axis.

    Keyword Args:
        axes:           :py:class:`~matplotlib.axes.Axes` object. If :py:obj:`None`,
          a new figure is created.
        xlabel:         x-axis label. Default: variable name.
        ylabel:         y-axis label. Default: observable axis label.

    Additional keyword arguments are transmitted to the :py:class:`~matplotlib.axes.Axes` creation function.

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
        >>> plot_response(var, obsleft, np.arange(3.0e9, 6.01e9, 0.5e9))
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
        ...             name=r"$\beta_x$",
        ...             plot_fmt={"linewidth": 3.0, "marker": "o"}
        ...         ),
        ...         at.LocalOpticsObservable(
        ...             [0], "beta", plane="y", name=r"$\beta_z$", plot_fmt="--"
        ...         )
        ...     ],
        ...     ring=ring
        ... )
        >>>
        >>> obsright =at.ObservableList(
        ...     [
        ...         at.GlobalOpticsObservable("tune", plane="x", name=r"$\nu_x$"),
        ...         at.GlobalOpticsObservable("tune", plane="y", name=r"$\nu_x$"),
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

        Example varying an evaluation parameter:

        >>> obs =at.ObservableList(
        ...     [
        ...         at.LocalOpticsObservable([0], "beta", plane="x", name=r"$\beta_x$"),
        ...         at.LocalOpticsObservable([0], "beta", plane="y", name=r"$\beta_z$")
        ...     ],
        ...     ring=ring,
        ...     dp = 0.0
        ... )
        >>> var = at.EvaluationVariable(obsleft, "dp", name=r"$\delta$")
        >>> ax=at.plot_response(
        ... var, obsleft, np.arange(-0.03, 0.0301,0.001), ylabel=r"$\beta\;[m]$"
        ... )

        .. image:: /images/delta_response.*
           :alt: delta response
    """

    def compute(v, obs):
        """Evaluate the observables for 1 variable value."""
        var.value = v
        for ob in obs:
            yield from ob.evaluate()

    def axes1(axes: Axes, obs: ObservableList, ylabel: str):
        """Plot all observables on a given axis."""

        def plot1(obs, ncurve):
            """Plot 1 curve."""
            fmt = getattr(obs, "plot_fmt", f"C{ncurve}")
            if isinstance(fmt, Mapping):
                return axes.plot(xx, next(values), label=obs.label, **fmt)
            else:
                return axes.plot(xx, next(values), fmt, label=obs.label)

        axes.set_ylabel(ylabel or obs.axis_label)
        for ob in obs:
            yield from plot1(ob, next(line_counter))

    if isinstance(args[0], ObservableList):
        obsleft, rng, *obsright = args
    else:
        rng, obsleft, *obsright = args

    if axes is None:
        _, axleft = plt.subplots(subplot_kw=kwargs)
    elif isinstance(axes, Axes):
        axleft = axes
    else:
        msg = "The 'axes' argument must be an Axes object."
        raise ValueError(msg)

    allaxes = (axleft, axleft.twinx()) if obsright else (axleft,)
    allobs = (obsleft, *obsright)

    # Compute all the observable values
    with (var.restore()):
        vals = [(v, *compute(v, allobs)) for v in rng]

    xx, *yy = zip(*vals, strict=True)
    values = iter(yy)
    line_counter = itertools.count()

    lines = []
    for ax, obs, ylab in zip(allaxes, allobs, (ylabel, ""), strict=False):
        lines.extend(axes1(ax, obs, ylab))

    axleft.set_xlabel(xlabel or var.name)
    axleft.legend(handles=lines)
    axleft.grid(True)

    return allaxes
