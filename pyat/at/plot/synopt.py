"""Lattice synoptics."""

from __future__ import annotations

__all__ = ["plot_synopt"]

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

from ..lattice import Lattice, Refpts, elements as elts

# Default properties for element representation
DIPOLE = {"label": "Dipoles", "facecolor": (0.5, 0.5, 1.0)}
QUADRUPOLE = {"label": "Quadrupoles", "facecolor": (1.0, 0.5, 0.5)}
SEXTUPOLE = {"label": "Sextupoles", "facecolor": (0.5, 1.0, 0.5)}
MULTIPOLE = {"label": "Multipoles", "facecolor": (0.25, 0.75, 0.25)}
MONITOR = {"label": "Monitors", "linestyle": None, "marker": 10, "color": "k"}


# noinspection PyDefaultArgument
def plot_synopt(
    ring: Lattice,
    axes: matplotlib.axes.Axes = None,
    dipole: dict | None = {},  # noqa: B006
    quadrupole: dict | None = {},  # noqa: B006
    sextupole: dict | None = {},  # noqa: B006
    multipole: dict | None = {},  # noqa: B006
    monitor: dict | None = {},  # noqa: B006
    labels: Refpts = None,
    **kwargs,
):
    """Plots a synoptic of a lattice.

    Parameters:
        ring:           Lattice description.
        axes:           :py:class:`~matplotlib.axes.Axes` for plotting the
          synoptic. If :py:obj:`None`, a new figure will be created. Otherwise,
          a new axes object sharing the same x-axis as the given one is created.
        labels (Refpts):    Select the elements for which the name is displayed.
          Default: :py:obj:`None`,
        dipole (dict):      Dictionary of properties overloading the default
          properties of the dipole representation.
          Example: :pycode:`{"facecolor": "xkcd:electric blue"}`. If :py:obj:`None`,
          dipoles are not shown.
        quadrupole (dict):  Same definition as for dipole,
        sextupole (dict):   Same definition as for dipole,
        multipole (dict):   Same definition as for dipole,
        monitor (dict):     Same definition as for dipole.

    Returns:
        synopt_axes (Axes): Synoptic axes
    """

    class Dipole(Polygon):
        xx = np.array([0, 0, 1, 1], dtype=float)
        yy = np.array([0, 1, 1, 0], dtype=float)

        def __init__(self, s, length, **kwargs):
            xy = np.stack((self.xx * length + s, self.yy), axis=1)
            super().__init__(xy, closed=False, **kwargs)

    class Quadrupole(Polygon):
        xx = np.array([0, 0, 0.5, 1, 1])
        yy = {
            True: np.array([0, 1, 1.4, 1, 0]),
            False: np.array([0, 1, 0.6, 1, 0]),
        }

        def __init__(self, s, length, foc, **kwargs):
            xy = np.stack((self.xx * length + s, self.yy[foc]), axis=1)
            super().__init__(xy, closed=False, **kwargs)

    class Sextupole(Polygon):
        xx = np.array([0, 0, 0.33, 0.66, 1, 1])
        yy = {
            True: np.array([0, 0.8, 1, 1, 0.8, 0]),
            False: np.array([0, 0.8, 0.6, 0.6, 0.8, 0]),
        }

        def __init__(self, s, length, foc, **kwargs):
            xy = np.stack((self.xx * length + s, self.yy[foc]), axis=1)
            super().__init__(xy, closed=False, **kwargs)

    class Multipole(Polygon):
        xx = np.array([0, 0, 1, 1], dtype=float)
        yy = np.array([0, 0.8, 0.8, 0])

        def __init__(self, s, length, **kwargs):
            xy = np.stack((self.xx * length + s, self.yy), axis=1)
            super().__init__(xy, closed=False, **kwargs)

    class Monitor(Polygon):
        xx = np.array([0.0, 0.0])
        yy = np.array([0.0, 1.2])

        def __init__(self, s, **kwargs):
            xy = np.stack((self.xx + s, self.yy), axis=1)
            super().__init__(xy, closed=False, **kwargs)

    def ismultipole(elem):
        return isinstance(elem, elts.Multipole) and not isinstance(
            elem, (elts.Dipole, elts.Quadrupole, elts.Sextupole)
        )

    if axes is None:
        _, axes = plt.subplots(subplot_kw=kwargs)

    axes.set_facecolor((0.0, 0.0, 0.0, 0.0))
    axsyn = axes.twinx()
    axsyn.set_title(ring.name, fontdict={"fontsize": "medium"}, loc="left")
    axsyn.set_axis_off()  # Set axis invisible
    axsyn.set_xlim(ring.s_range)
    axsyn.set_ylim((0.0, 20.0))  # Initial scaling of elements
    axsyn.set_zorder(-0.2)  # Put synoptic in the background

    s_pos = ring.get_s_pos(range(len(ring)))

    if dipole is not None:
        props = DIPOLE | dipole
        dipoles = PatchCollection(
            (
                Dipole(s, el.Length)
                for s, el in zip(s_pos, ring, strict=True)
                if isinstance(el, elts.Dipole)
            ),
            **props,
        )
        axsyn.add_collection(dipoles)

    if quadrupole is not None:
        props = QUADRUPOLE | quadrupole
        quadrupoles = PatchCollection(
            (
                Quadrupole(s, el.Length, el.PolynomB[1] >= 0.0)
                for s, el in zip(s_pos, ring, strict=True)
                if isinstance(el, elts.Quadrupole)
            ),
            **props,
        )
        axsyn.add_collection(quadrupoles)

    if sextupole is not None:
        props = SEXTUPOLE | sextupole
        sextupoles = PatchCollection(
            (
                Sextupole(s, el.Length, el.PolynomB[2] >= 0.0)
                for s, el in zip(s_pos, ring, strict=True)
                if isinstance(el, elts.Sextupole)
            ),
            **props,
        )
        axsyn.add_collection(sextupoles)

    if multipole is not None:
        props = MULTIPOLE | multipole
        multipoles = PatchCollection(
            (
                Multipole(s, el.Length)
                for s, el in zip(s_pos, ring, strict=True)
                if ismultipole(el)
            ),
            **props,
        )
        axsyn.add_collection(multipoles)

    if monitor is not None:
        props = MONITOR | monitor
        s = s_pos[[isinstance(el, elts.Monitor) for el in ring]]
        y = np.zeros(s.shape)
        axsyn.plot(s, y, **props)

    for idx in ring.get_uint32_index(labels):
        el = ring[idx]
        s = s_pos[idx]
        axsyn.text(
            s + 0.5 * el.Length,
            1.6,
            el.FamName[:10],
            rotation="vertical",
            horizontalalignment="center",
            fontsize="small",
        )

    return axsyn, axes
