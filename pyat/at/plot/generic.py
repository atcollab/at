"""AT generic plotting function"""

from __future__ import annotations

__all__ = ["baseplot"]

from itertools import chain, repeat
from collections.abc import Callable

# noinspection PyPackageRequirements
import matplotlib.pyplot as plt

from .synopt import plot_synopt
from ..lattice import Lattice

SLICES = 400


def baseplot(ring: Lattice, plot_function: Callable, *args, **kwargs):
    """Generic lattice plot

    :py:func:`baseplot` divides the region of interest of ring into small
    elements, calls the specified function to get the plot data and calls
    matplotlib functions to generate the plot.
    By default, it creates a new figure for the plot, but if provided with
    :py:class:`~matplotlib.axes.Axes` objects it can be used as part of a GUI

    Parameters:
        ring:           Lattice description.
        plot_function:  Specific data generating function to be called
          plotting function. ``plot_function`` is called as:

          :code:`title, left, right = plot_function(ring, refpts, *args,
          **kwargs)`

          and should return 2 or 3 output:

          ``title``: plot title or :py:obj:`None`

          ``left``: tuple returning the data for the main (left) axis

            left[0] - y-axis label

            left[1] - xdata: (N,) array (s coordinate)

            left[2] - ydata: iterable of (N,) or (N,M) arrays. Lines from a
            (N, M) array share the same style and label

            left[3]   labels: (optional) iterable of strings as long as ydata

          ``right``: tuple returning the data for the secondary (right) axis
        *args:          All other positional parameters are sent to the plotting
          function

    Keyword Args:
        s_range:            Lattice range of interest, default: unchanged,
          initially set to the full cell.
        axes (tuple[Axes, Optional[Axes]): :py:class:`~matplotlib.axes.Axes`
          for plotting as (primary_axes, secondary_axes).
          Default: create new axes
        slices (int):       Number of slices. Default: 400
        legend (bool):      Show a legend on the plot. Default: :py:obj:`True`
        labels (Refpts):    display the name of selected elements.
          Default: :py:obj:`None`
        block (bool):       If :py:obj:`True`, block until the figure is closed.
          Default: :py:obj:`False`
        dipole (dict):      Dictionary of properties overloading the default
          properties of dipole representation. See :py:func:`.plot_synopt`
          for details
        quadrupole (dict):  Same definition as for dipole
        sextupole (dict):   Same definition as for dipole
        multipole (dict):   Same definition as for dipole
        monitor (dict):     Same definition as for dipole
        **kwargs:           All other keywords are sent to the plotting function

    Returns:
        left_axes (Axes):   Main (left) axes
        right_axes (Axes):  Secondary (right) axes or :py:obj:`None`
        synopt_axes (Axes): Synoptic axes
    """

    def get_synopt(ax, **kw):
        axs = getattr(ax, "axsyn", None)
        if axs is None:
            axs = plot_synopt(ring, axes=ax, **kw)
            ax.axsyn = axs
        return axs

    def get_axright(ax):
        axr = getattr(ax, "axright", None)
        if axr is None:
            axr = axleft.twinx()
            ax.axright = axr
            # axright._get_lines = ax._get_lines
        return axr

    def get_style(ax):
        style_iter = getattr(ax, "style", None)
        if style_iter is None:
            style_iter = iter(plt.rcParams["axes.prop_cycle"])
            ax.style = style_iter
        return style_iter

    def plot1(ax, props, yaxis_label, x, y, labels=()):
        # for y1, label in zip(y, chain(labels, repeat(None))):
        #     ax.plot(x, y1, label=label)
        for y1, prop, label in zip(y, props, chain(labels, repeat(None))):
            ax.plot(x, y1, label=label, **prop)
        ax.set_ylabel(yaxis_label)

    def new_axes(ax, data):
        ax.set(
            xlim=ring.s_range,
            xlabel="s [m]",
            facecolor=[1.0, 1.0, 1.0, 0.0],
            title=title,
        )
        ax.set_title(ring.name, fontdict={"fontsize": "medium"}, loc="left")
        props = get_style(ax)
        axs = get_synopt(ax, **synargs)
        plot1(ax, props, *data)
        ax.grid(visible=True)
        return props, axs

    # extract baseplot arguments
    slices = kwargs.pop("slices", SLICES)
    axes = kwargs.pop("axes", None)
    legend = kwargs.pop("legend", True)
    block = kwargs.pop("block", False)
    s_range = kwargs.pop("s_range", None)

    if s_range is not None:
        ring.s_range = s_range

    try:
        axleft, axright = axes
    except TypeError:
        axleft = axes
        axright = None

    # extract synopt arguments
    synkeys = ["dipole", "quadrupole", "sextupole", "multipole", "monitor", "labels"]
    kwkeys = list(kwargs.keys())
    synargs = {k: kwargs.pop(k) for k in kwkeys if k in synkeys}

    # slice the ring
    rg = ring.slice(slices=slices)

    # get the data for the plot
    title, *plots = plot_function(rg, rg.i_range, *args, **kwargs)

    # prepare the axes
    if axleft is None:
        # Create new axes
        fig, axleft = plt.subplots()

    lprops, axsyn = new_axes(axleft, plots[0])
    handles, labs = axleft.get_legend_handles_labels()

    if len(plots) > 1:
        if axright is None:
            # Add a right axis on axleft
            axright = get_axright(axleft)
            plot1(axright, lprops, *plots[1])
            h2, l2 = axright.get_legend_handles_labels()
            handles.extend(h2)
            labs.extend(l2)
        elif axleft.get_shared_x_axes().joined(axleft, axright):
            # Already right axis of axleft
            plot1(axright, lprops, *plots[1])
            h2, l2 = axright.get_legend_handles_labels()
            handles.extend(h2)
            labs.extend(l2)
        else:
            # New axes
            _, _ = new_axes(axright, plots[1])
            if legend:
                axright.legend()

    if legend and handles:
        axleft.legend(handles, labs)

    plt.show(block=block)
    return axleft, axright, axsyn
