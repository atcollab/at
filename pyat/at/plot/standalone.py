"""AT plotting functions"""
from typing import Tuple
from at.lattice import Lattice
import matplotlib.pyplot as plt
import numpy
from numpy import ndarray


# Function to compute and plot acceptance
def plot_acceptance(ring: Lattice, *args, **kwargs):
    # noinspection PyUnresolvedReferences
    r"""Plots the acceptance

    Computes the acceptance at repfts observation points using
    :py:func:`.get_acceptance` and plots the tracked
    and survived particles, and the acceptance boundary.

    Parameters:
        ring:           Lattice definition

    Keyword Args:
        acceptance(Tuple[ndarray, ndarray, ndarray]): Tuple containing
          pre-computed acceptance ``(boundary, survived, grid)``
        planes:         max. dimension 2, Plane(s) to scan for the acceptance.
          Allowed values are: ``'x'``, ``'xp'``, ``'y'``,
          ``'yp'``, ``'dp'``, ``'ct'``
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
          enabled, ``grid_mode`` is forced to :py:attr:`.GridMode.RECURSIVE`
          (most efficient in single core)
        verbose (bool):     Print out some information
        divider (int):      Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero:
        start_method (str): Python multiprocessing start method. The default
          :py:obj:`None` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for macOS and
          Windows is ``'spawn'``. ``'fork'`` may used for macOS to speed up
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

    units = {'x': '[m]', 'xp': '[rad]', 'y': '[m]',
             'yp': '[rad]', 'dp': '', 'ct': '[m]'}

    obspt = kwargs.pop('obspt', None)
    block = kwargs.pop('block', False)
    acceptance = kwargs.pop('acceptance', None)
    if obspt is not None:
        assert numpy.isscalar(obspt), 'Scalar value needed for obspt'
    kwargs['refpts'] = obspt
    if acceptance is None:
        boundary, survived, grid = ring.get_acceptance(*args, **kwargs)
    else:
        boundary, survived, grid = acceptance
    planes = args[0]
    plt.figure()
    plt.plot(*grid, '.', label='Tracked particles')
    plt.plot(*survived, '.', label='Survived particles')
    if len(planes) == 1:
        plt.plot(boundary, numpy.zeros(2), label='Acceptance')
        plt.title('1D {0} acceptance'.format(planes[0]))
        plt.xlabel('{0} {1}'.format(planes[0], units[planes[0]]))
    else:
        plt.plot(*boundary, label='Acceptance')
        plt.title('2D {0} {1} acceptance'.format(planes[0], planes[1]))
        plt.xlabel('{0} {1}'.format(planes[0], units[planes[0]]))
        plt.ylabel('{0} {1}'.format(planes[1], units[planes[1]]))
    plt.legend()
    plt.show(block=block)
    return boundary, survived, grid


def plot_geometry(ring: Lattice, start_coordinates=(0, 0, 0),
                  centered=False, offset=(0.0, 0.0), ax=None, label=''):
    """Compute the 2D ring geometry in cartesian coordinates

    Parameters:
        ring: Lattice description
        start_coordinates: x,y,angle at starting point
        centered: it True the coordinates origin is the center of the ring
        offset: (dx, dy) offsets coordinates by the given amount
        ax: axes where to plot, if not given axes are created
        label: label of curve

    Returns:
        geomdata: recarray containing, x, y, angle
        radius: machine radius
        ax: plot axis
    """
    if not ax:
        fig, ax = plt.subplots()
    geom, radius = ring.get_geometry(start_coordinates=start_coordinates,
                                     centered=centered, offset=offset)
    ax.plot(geom['x'], geom['y'], 'o:', linewidth=0.5, markersize=2,
            label=label)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_aspect('equal', 'box')
    return geom, radius, ax


Lattice.plot_acceptance = plot_acceptance
Lattice.plot_geometry = plot_geometry
