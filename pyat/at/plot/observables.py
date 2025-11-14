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


class _RingSplitter:
    """generate a ring split in small slices for plotting"""

    def __init__(self, ring: Lattice, s_range: tuple[float, float], slices: int):
        ring.s_range = s_range
        self._ring0 = ring
        self._ring1 = None
        self.s = None
        self.slices = slices

    @property
    def ring(self):
        if self._ring1:
            return self._ring1
        else:
            rg = self._ring0.slice(slices=self.slices)
            self.s = rg.get_s_pos(range(len(rg) + 1))
            self._ring1 = rg
            return rg


def plot_observables(
        ring: Lattice,
        obsleft: ObservableList,
        *obsright: ObservableList,
        s_range: tuple[float, float] | None = None,
        axes: Axes | None = None,
        slices: int = _SLICES,
        **kwargs,
) -> tuple[Axes]:
    # noinspection PyUnresolvedReferences
    r"""Plot element observables along a lattice.

    Args:
        ring:   Lattice description
        obsleft: List of :py:class:`.ElementObservable` plotted against the left axis.
          if refpts is :py:obj:`.All`, a line is drawn. Otherwise, markers are drawn.
          It is recommended to use Observables with scalar values. Otherwise, all the
          values are plotted but share the same line properties and legend,
        obsright: Optional list of :py:class:`.ElementObservable` plotted against the
          right axis,
        axes: :py:class:`~.matplotlib.axes.Axes` in which the observables are plotted.
          if :py:obj:`None`, a new figure is created,
        s_range:            Lattice range of interest, default: whole lattice,
        slices: Number of lattice slices for getting smooth curves. Default: 400.

    The following keywords are transmitted to the :py:func:`.plot_synopt` function.They
    apply to the main (left) axis and are ignored when plotting in exising axes:

    Keyword Args:
        labels (Refpts):    Select elements for which the name is displayed.
          Default: :py:obj:`None`,
        dipole (dict):      Dictionary of properties overloading the default
          properties of the dipole representation.
          Example: :pycode:`{"facecolor": "xkcd:electric blue"}`. If :py:obj:`None`,
          dipoles are not shown.
        quadrupole (dict):  Same definition as for dipole,
        sextupole (dict):   Same definition as for dipole,
        multipole (dict):   Same definition as for dipole,
        monitor (dict):     Same definition as for dipole.

    The following keyword arguments are transmitted to the
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

    Examples:
        Minimal example using default values:

        >>> obsmu = at.ObservableList(
        ...     [
        ...         at.LocalOpticsObservable(at.All, "mu", plane=0),
        ...         at.LocalOpticsObservable(at.All, "mu", plane=1),
        ...     ]
        ... )
        >>>
        >>> ax1, = at.plot_observables(ring, obsmu)

        .. image:: /images/plot_observables/phase_obs.*
           :alt: phase advance plot

        This example demonstrates the use of the postfun post-processing attribute of
        observables to plot the beam envelopes for arbitrary emitttance values:

        >>> # Define the emittances
        >>> emit_x = 130.0e-12
        >>> emit_y = 10.0e-12
        >>>
        >>> # beam size computation
        >>> sigma_x = lambda x: 1.0e6 * np.sqrt(x * emit_x)  # result in um
        >>> sigma_y = lambda y: 1.0e6 * np.sqrt(y * emit_y)  # result in um
        >>>
        >>> # Observables
        >>> obsenv = at.ObservableList(
        ...     [
        ...         at.LocalOpticsObservable(
        ...             at.All, "beta", plane="x", postfun=sigma_x, label=r"$\sigma_x$"
        ...         ),
        ...         at.LocalOpticsObservable(
        ...             at.All, "beta", plane="y", postfun=sigma_y, label=r"$\sigma_y$"
        ...         ),
        ...     ]
        ... )
        >>> (ax2,) = at.plot_observables(ring, obsenv, title="Beam envelopes")
        >>> ax1.set_ylabel(r"Beam size [${\mu}m$]")

        .. image:: /images/plot_observables/envelope_obs.*
           :alt: envelope plot
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
                x = curve_s
            else:
                args, kwargs = get_format(obs, f"oC{next(dot_counter)}")
                x = dot_s[obs._boolrefs]

            return axes.plot(x, obs.value, *args, label=obs.label, **kwargs)

        axes.set_ylabel(obs.axis_label)
        for ob in obs:
            yield from plot1(ob)

    def evaluate(obs: ObservableList) -> None:
        """Evaluates one observable."""
        curves = ObservableList(**obs.kwargs)
        dots = ObservableList(**obs.kwargs)

        for ob in obs:
            if not isinstance(ob, ElementObservable):
                msg = f"{ob.name} is not an ElementObservable object."
                raise ValueError(msg)
            elif ob.refpts is All:
                curves.append(ob)
            else:
                dots.append(ob)

        if curves:
            # Evaluate curve data
            curves.evaluate(ring=splitter.ring)

        if dots:
            # Evaluate marker data
            dots.evaluate(ring=ring)

    splitter = _RingSplitter(ring, s_range, slices)

    if axes is None:
        _, axleft = plot_synopt(ring, **kwargs)
    elif isinstance(axes, Axes):
        axleft = axes
    else:
        msg = "The 'axes' argument must be an Axes object."
        raise ValueError(msg)

    allaxes = (axleft, axleft.twinx()) if obsright else (axleft,)
    allobs = (obsleft, *obsright)

    for obs in allobs:
        evaluate(obs)

    curve_s = splitter.s
    dot_s = ring.get_s_pos(range(len(ring) + 1))
    curve_counter = itertools.count()
    dot_counter = itertools.count()

    lines = []
    for ax, obs in zip(allaxes, allobs, strict=True):
        lines.extend(axes1(ax, obs))

    axleft.set_xlabel("s [m]")
    axleft.legend(handles=[ln for ln in lines if not ln.get_label().startswith("_")])
    axleft.grid(True)
    return allaxes
