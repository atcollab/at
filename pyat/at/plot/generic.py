"""AT generic plotting function"""
from __future__ import annotations
from itertools import chain, repeat
# noinspection PyPackageRequirements
import matplotlib.pyplot as plt
from typing import Callable
from .synopt import plot_synopt
from ..lattice import Lattice

SLICES = 400

__all__ = ['baseplot']


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
        *args:          All other positional parameters are sent to the

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

    def plot1(ax, yaxis_label, x, y, labels=()):
        lines = []
        for y1, prop, label in zip(y, props, chain(labels, repeat(None))):
            ll = ax.plot(x, y1, **prop)
            if label is not None:
                ll[0].set_label(label)
            lines += ll
        ax.set_ylabel(yaxis_label)
        return lines

    def labeled(line):
        return not line.properties()['label'].startswith('_')

    # extract baseplot arguments
    slices = kwargs.pop('slices', SLICES)
    axes = kwargs.pop('axes', None)
    legend = kwargs.pop('legend', True)
    block = kwargs.pop('block', False)
    if 's_range' in kwargs:
        ring.s_range = kwargs.pop('s_range')

    # extract synopt arguments
    synkeys = ['dipole', 'quadrupole', 'sextupole', 'multipole',
               'monitor', 'labels']
    kwkeys = list(kwargs.keys())
    synargs = dict((k, kwargs.pop(k)) for k in kwkeys if k in synkeys)

    # get color cycle
    cycle_props = plt.rcParams['axes.prop_cycle']

    # slice the ring
    rg = ring.slice(slices=slices)

    # get the data for the plot
    pout = plot_function(rg, rg.i_range, *args, **kwargs)
    title = pout[0]
    plots = pout[1:]

    # prepare the axes
    if axes is None:
        # Create new axes
        nplots = len(plots)
        fig = plt.figure()
        axleft = fig.add_subplot(111, xlim=rg.s_range, xlabel='s [m]',
                                 facecolor=[1.0, 1.0, 1.0, 0.0],
                                 title=title)
        axright = axleft.twinx() if (nplots >= 2) else None
        axleft.set_title(ring.name, fontdict={'fontsize': 'medium'},
                         loc='left')
        axsyn = plot_synopt(ring, axes=axleft, **synargs)
    else:
        # Use existing axes
        axleft, axright = axes
        axsyn = None
        nplots = 1 if axright is None else len(plots)

    props = iter(cycle_props())

    # left plot
    lines1 = plot1(axleft, *plots[0])
    # right plot
    lines2 = [] if (nplots < 2) else plot1(axright, *plots[1])
    if legend:
        if nplots < 2:
            axleft.legend(handles=[li for li in lines1 if labeled(li)])
        elif axleft.get_shared_x_axes().joined(axleft, axright):
            axleft.legend(
                handles=[li for li in lines1 + lines2 if labeled(li)])
        else:
            axleft.legend(handles=[li for li in lines1 if labeled(li)])
            axright.legend(handles=[li for li in lines2 if labeled(li)])
    plt.show(block=block)
    return axleft, axright, axsyn
