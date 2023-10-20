"""AT plotting functions"""
from __future__ import annotations
from at.lattice import Lattice, axis_
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy
from numpy import ndarray
from math import sqrt


# Function to compute and plot acceptance
def plot_acceptance(ring: Lattice, planes, *args, **kwargs):
    # noinspection PyUnresolvedReferences
    r"""Plots the acceptance

    Computes the acceptance at repfts observation points using
    :py:func:`.get_acceptance` and plots the tracked
    and survived particles, and the acceptance boundary.

    Parameters:
        ring:           Lattice definition
        planes:         max. dimension 2, Plane(s) to scan for the acceptance.
          Allowed values are: *'x'*, *'xp'*, *'y'*,
          *'yp'*, *'dp'*, *'ct'*

    Keyword Args:
        acceptance (tuple[ndarray, ndarray, ndarray]): tuple containing
          pre-computed acceptance *(boundary, survived, grid)*
        npoints:        (len(planes),) array: number of points in each
          dimension
        amplitudes:     (len(planes),) array: set the search range:

          * :py:attr:`GridMode.CARTESIAN/RADIAL <.GridMode.RADIAL>`:
            max. amplitude
          * :py:attr:`.GridMode.RECURSIVE`: initial step
        nturns (int):       Number of turns for the tracking
        obspt (Refpts):    Observation points. Default: start of the machine
        dp (float):         Static momentum offset
        offset:             Initial orbit. Default: closed orbit
        bounds:             Defines the tracked range: range=bounds*amplitude.
          It can be used to select quadrants. For example, default values are:

          * :py:attr:`.GridMode.CARTESIAN`: ((-1, 1), (0, 1))
          * :py:attr:`GridMode.RADIAL/RECURSIVE <.GridMode.RADIAL>`: ((0, 1),
            (:math:`\pi`, 0))
        grid_mode (GridMode):   Defines the evaluation grid:

          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp (bool):      Use python multiprocessing (:py:func:`.patpass`,
          default use :py:func:`.lattice_pass`). In case multiprocessing is not
          enabled, *grid_mode* is forced to :py:attr:`.GridMode.RECURSIVE`
          (most efficient in single core)
        verbose (bool):     Print out some information
        divider (int):      Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero:
        start_method (str): Python multiprocessing start method. The default
          :py:obj:`None` uses the python default that is considered safe.
          Available parameters: *'fork'*, *'spawn'*, *'forkserver'*.
          The default for linux is *'fork'*, the default for macOS and
          Windows is *'spawn'*. *'fork'* may be used for macOS to speed up
          the calculation or to solve runtime errors, however  it is
          considered unsafe.

    Returns:
        boundary:   (2,n) array: 2D acceptance
        tracked:    (2,n) array: Coordinates of tracked particles
        survived:   (2,n) array: Coordinates of surviving particles

    Example:
        >>> ring.plot_acceptance(planes, npoints, amplitudes)
        >>> plt.show()
    """
    obspt = kwargs.pop('obspt', None)
    block = kwargs.pop('block', False)
    acceptance = kwargs.pop('acceptance', None)
    if obspt is not None:
        assert numpy.isscalar(obspt), 'Scalar value needed for obspt'
    kwargs['refpts'] = obspt
    if acceptance is None:
        boundary, survived, grid = ring.get_acceptance(planes, *args, **kwargs)
    else:
        boundary, survived, grid = acceptance
    plt.figure()
    plt.plot(*grid, '.', label='Tracked particles')
    plt.plot(*survived, '.', label='Survived particles')
    if len(planes) == 1:
        pl0 = axis_(planes[0])
        plt.plot(boundary, numpy.zeros(2), label='Acceptance')
        plt.title('1D {0} acceptance'.format(pl0['label']))
        plt.xlabel('{0}{1}'.format(pl0['label'], pl0['unit']))
    else:
        pl0, pl1 = axis_(planes)
        plt.plot(*boundary, label='Acceptance')
        plt.title('2D {0}-{1} acceptance'.format(pl0['label'], pl1['label']))
        plt.xlabel('{0}{1}'.format(pl0['label'], pl0['unit']))
        plt.xlabel('{0}{1}'.format(pl1['label'], pl1['unit']))
    plt.legend()
    plt.show(block=block)
    return boundary, survived, grid


def plot_geometry(ring: Lattice,
                  start_coordinates: tuple[float, float, float] = (0, 0, 0),
                  centered: bool = False, ax: Axes = None, **kwargs):
    """Compute and plot the 2D ring geometry in cartesian coordinates.

    Parameters:
        ring: Lattice description
        start_coordinates: x,y,angle at starting point
        centered: it True the coordinates origin is the center of the ring
        ax: axes for the plot. If not given, new axes are created

    Keyword arguments are forwarded to the underlying
    :py:func:`~matplotlib.pyplot.plot` function

    Returns:
        geomdata: recarray containing, x, y, angle
        radius: machine radius
        ax: plot axis

    Example:
        >>> ring.plot_geometry()
    """
    if not ax:
        fig, ax = plt.subplots()
    geom, radius = ring.get_geometry(start_coordinates=start_coordinates,
                                     centered=centered)
    ax.plot(geom['x'], geom['y'], 'o:',
            linewidth=kwargs.pop('linewidth', 0.5),
            markersize=kwargs.pop('markersize', 2),
            **kwargs)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_aspect('equal', 'box')
    return geom, radius, ax


def plot_sigma(sigma, axis: tuple[str, str] = ('x', 'xp'), scale: float = 1.0,
               ax: Axes = None, **kwargs):
    r"""Plot the projection of the phase space defined by a
    :math:`\Sigma`-matrix on the selected plane.

    Arguments:
        sigma:  :math:`\Sigma`-matrix
        axis:   tuple if indices defining the plane of the :math:`\Sigma`
          projection. Allowed values are: *'x'*, *'xp'*, *'y'*,
          *'yp'*, *'dp'*, *'ct'*. Default: (*'x'*, *'xp'*)
        scale:  Scaling factor for the emittance
        ax: axes for the plot. If not given, new axes are created

    Keyword arguments are forwarded to the underlying
    :py:func:`~matplotlib.pyplot.plot` function
    """
    if not ax:
        fig, ax = plt.subplots()
    ax1, ax2 = axis_(axis)
    axid = axis_(axis, key='index')
    sig22 = sigma[numpy.ix_(axid, axid)]
    eps = sqrt(sig22[0, 0] * sig22[1, 1] - sig22[1, 0] * sig22[0, 1])
    sigx = sqrt(sig22[0, 0])
    tr = numpy.array([[sigx, 0.0],
                      [sig22[0, 1] / sigx, eps / sigx]])
    loop = 2.0 * numpy.pi * numpy.arange(0.0, 1.0, 0.001)
    normcoord = numpy.vstack((numpy.cos(loop), numpy.sin(loop)))
    coord = tr @ normcoord
    line = ax.plot(scale*coord[0, :], scale*coord[1, :], **kwargs)
    ax.set_title('{0}-{1} phase space'.format(ax1['label'], ax2['label']))
    ax.set_xlabel('{0}{1}'.format(ax1['label'], ax1['unit']))
    ax.set_ylabel('{0}{1}'.format(ax2['label'], ax2['unit']))
    return line


Lattice.plot_acceptance = plot_acceptance
Lattice.plot_geometry = plot_geometry
