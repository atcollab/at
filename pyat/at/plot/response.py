"""Plot :py:class:`.Observable` values as a function of a
:py:class:`Variable <.VariableBase>`.
"""

from __future__ import annotations

__all__ = ["plot_response"]

from collections.abc import Mapping
import itertools
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from ..lattice import VariableBase
from ..latticetools import ObservableList


def plot_response(
    var: VariableBase,
    rng: Iterable[float],
    obsleft: ObservableList,
    *obsright: ObservableList,
    axes: Axes | None = None,
    xlabel: str = "",
    ylabel: str = "",
    **kwargs
) -> tuple[Axes]:
    # noinspection PyUnresolvedReferences
    r"""Plot :py:class:`.Observable` values as a function of a
    :py:class:`Variable <.VariableBase>`.

    Args:
        var:        Variable object,
        rng:        range of variation for the variable,
        obsleft:    List of Observables plotted on the left axis. It is recommended to
          use Observables with scalar values. Otherwise, all the values are plotted but
          share the same line properties and legend,
        obsright:   Optional list of Observables plotted on the right axis.
        axes:           :py:class:`~matplotlib.axes.Axes` object in which the figure
          is plotted. If :py:obj:`None`, a new figure is created.
        xlabel:         x-axis label. May contain Latex math code.
          Default: variable name.
        ylabel:         y-axis label. May contain Latex math code.
          Default: observable :py:attr:`~.ObservableList.axis_label`.

    Additional keyword arguments are transmitted to the
    :py:class:`~matplotlib.axes.Axes` creation function.They apply to the main (left)
    axis and are ignored when plotting in exising axes:

    Keyword Args:
        title (str):    Plot title,
        ylim (tuple):     Y-axis limits,
        *: for other keywords see
          :py:obj:`~.matplotlib.figure.Figure.add_subplot`


    Returns:
        axes: tuple of :py:class:`~.matplotlib.axes.Axes`. Contains 2 elements if there
          is a plot on the right y-axis, 1 element otherwise.

    Example:
        Minimal example using only default values:

        >>> obsl = at.ObservableList(
        ...     [at.EmittanceObservable("emittances", plane="x")],
        ...     ring=ring,
        ... )
        >>> obsr = at.ObservableList(
        ...     [at.EmittanceObservable("sigma_e")],
        ...     ring=ring,
        ... )
        >>> var = at.AttributeVariable(ring, "energy", name="energy [eV]")
        >>> ax1, ax2 = at.plot_response(
        ...     var, np.arange(3.0e9, 6.01e9, 0.5e9), obsl, obsr
        ... )
        >>>

        .. image:: /images/emittance_response.*
           :alt: emittance response

        Example showing the formatting possibilities by:

        - using the :py:attr:`.Observable.plot_fmt` attribute for line formatting,
        - using dual y-axis,
        - using the *ylim* and *title* parameters.

        >>> obsleft =at.ObservableList(
        ...     [
        ...         at.LocalOpticsObservable(
        ...             [0], "beta", plane="x",
        ...             plot_fmt={"linewidth": 3.0, "marker": "o"}
        ...         ),
        ...         at.LocalOpticsObservable([0], "beta", plane="y", plot_fmt="--")
        ...     ],
        ...     ring=ring
        ... )
        >>>
        >>> obsright =at.ObservableList(
        ...     [
        ...         at.GlobalOpticsObservable("tune", plane="x"),
        ...         at.GlobalOpticsObservable("tune", plane="y"),
        ...     ],
        ...     ring=ring
        ... )
        >>>
        >>> var = RefptsVariable(
        ...     "QF1[AE]", "Kn1L", name="QF1 integrated strength", ring=ring
        ... )
        >>> ax = at.plot_response(
        ...     var,
        ...     np.arange(0.732, 0.852, 0.01),
        ...     obsleft,
        ...     obsright,
        ...     ylim=[0.0, 10.0],
        ...     title="Example of plot_response"
        ... )

        .. image:: /images/beta_response.*
           :alt: beta response

        Example varying an evaluation parameter:

        >>> obs =at.ObservableList(
        ...     [
        ...         at.LocalOpticsObservable([0], "beta", plane="x"),
        ...         at.LocalOpticsObservable([0], "beta", plane="y")
        ...     ],
        ...     ring=ring,
        ...     dp = 0.0
        ... )
        >>> var = at.EvaluationVariable(obsleft, "dp", name=r"$\delta$")
        >>> ax=at.plot_response(var, np.arange(-0.03, 0.0301,0.001), obsleft)

        .. image:: /images/delta_response.*
           :alt: delta response
    """

    def compute(v, obs):
        """Evaluate the observables for 1 variable value."""
        var.value = v
        for ob in obs:
            yield from ob.evaluate()

    def axes1(axes: Axes, obs: ObservableList):
        """Plot all observables on a given axis."""

        def plot1(obs, ncurve):
            """Plot 1 curve."""
            fmt = getattr(obs, "plot_fmt", f"C{ncurve}")
            if isinstance(fmt, Mapping):
                return axes.plot(xx, next(values), label=obs.label, **fmt)
            else:
                return axes.plot(xx, next(values), fmt, label=obs.label)

        axes.set_ylabel(obs.axis_label)
        for ob in obs:
            yield from plot1(ob, next(line_counter))

    if isinstance(rng, ObservableList):
        # swap arguments for the old argument order
        rng, obsleft = obsleft, rng

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
    for ax, obs in zip(allaxes, allobs, strict=True):
        lines.extend(axes1(ax, obs))

    axleft.set_xlabel(xlabel or var.name)
    if ylabel:
        axleft.set_ylabel(ylabel)
    axleft.legend(handles=lines)
    axleft.grid(True)

    return allaxes
