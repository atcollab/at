from __future__ import annotations

__all__ = ["plot_observables"]

from collections.abc import Mapping, Generator
import itertools

from matplotlib.axes import Axes
from matplotlib.lines import Line2D

from ..lattice import Lattice, All
from ..latticetools import ElementObservable, ObservableList
from .synopt import plot_synopt

_SLICES = 400


def plot_observables(
    ring: Lattice,
    obsleft: ObservableList,
    *obsright: ObservableList,
    s_range: tuple[float, float] | None = None,
    axes: Axes | None = None,
    slices: int = _SLICES,
    **kwargs,
) -> tuple[Axes]:
    r"""Plot element observables along a lattice.

    Args:
        ring:   Lattice description
        obsleft: List of :py:class:`.ElementObservable` plotted against the left axis.
          if refpts is :py:obj:`.All`, a line is drawn. Otherwise, markers are drawn.
        obsright: List of :py:class:`.ElementObservable` plotted against the right axis.
        axes: :py:class:`~.matplotlib.axes.Axes` in which the observables are plotted.
          if :py:obj:`None`, new axes are created.
        s_range:            Lattice range of interest, default: whole lattice,
        slices: Number of lattice slices for getting smooth curves. Default: 400.

    Keyword Args:
        title (str):        Plot title,
        labels (Refpts):    Select elements for which the name is displayed.
          Default: :py:obj:`None`
        dipole (dict):      Dictionary of properties overloading the default
          properties of dipole representation. See :py:func:`.plot_synopt`
          for details
        quadrupole (dict):  Same definition as for dipole
        sextupole (dict):   Same definition as for dipole
        multipole (dict):   Same definition as for dipole
        monitor (dict):     Same definition as for dipole

    Returns:
        axes: tuple of :py:class:`~.matplotlib.axes.Axes`. Contains 1 element if no
          plot on the right axis, 2 elements otherwise.
    """

    def get_format(obs: ElementObservable, default: str) -> tuple[tuple, Mapping]:
        fmt = getattr(obs, "plot_fmt", default)
        if isinstance(fmt, Mapping):
            return (), fmt
        else:
            return (fmt,), {}

    def axes1(axes: Axes, obs: ObservableList) -> Generator[Line2D, None, None]:
        """Plot all observables on a given axis."""

        def plot1(obs: ElementObservable):
            """Plot a single observable."""
            if obs.refpts is All:
                args, kwargs = get_format(obs, f"C{next(curve_counter)}")
                return axes.plot(
                    curve_s, next(curve_values), *args, label=obs.label, **kwargs
                )
            else:
                args, kwargs = get_format(obs, f"oC{next(dot_counter)}")
                x = dot_s[obs._boolrefs]
                return axes.plot(x, next(dot_values), *args, label=obs.label, **kwargs)

        axes.set_ylabel(obs.axis_label)
        for ob in obs:
            yield from plot1(ob)

    def select(obs: ObservableList) -> None:
        """Distribute observables between curves and dots"""
        for ob in obs:
            if not isinstance(ob, ElementObservable):
                msg = f"{ob.name} is not an ElementObservable object."
                raise ValueError(msg)
            elif ob.refpts is All:
                curves.append(ob)
            else:
                dots.append(ob)
        dots.kwargs.update(obs.kwargs)
        curves.kwargs.update(obs.kwargs)

    ring.s_range = s_range

    if axes is None:
        _, axleft = plot_synopt(ring, **kwargs)
    elif isinstance(axes, Axes):
        axleft = axes
    else:
        msg = "The 'axes' argument must be an Axes object."
        raise ValueError(msg)

    allaxes = (axleft, axleft.twinx()) if obsright else (axleft,)
    allobs = (obsleft, *obsright)

    curves = ObservableList()
    dots = ObservableList()
    for obs in allobs:
        select(obs)

    if curves:
        # Evaluate curve data
        rg = ring.slice(slices=slices)
        curve_s = rg.get_s_pos(range(len(rg) + 1))
        curve_values = iter(curves.evaluate(ring=rg))

    if dots:
        # Evaluate marker data
        dot_s = ring.get_s_pos(range(len(ring) + 1))
        dot_values = iter(dots.evaluate(ring=ring))

    curve_counter = itertools.count()
    dot_counter = itertools.count()

    lines = []
    for ax, obs in zip(allaxes, allobs, strict=True):
        lines.extend(axes1(ax, obs))

    axleft.set_xlabel("s [m]")
    axleft.legend(handles=[l for l in lines if not l.get_label().startswith("_")])
    axleft.grid(True)
    return allaxes
