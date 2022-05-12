"""
Module to calculate the 2D or 1D acceptance of the ring
"""

import numpy
from .boundary import GridMode
from .boundary import boundary_search
import multiprocessing
from ..lattice import Lattice


__all__ = ['get_acceptance', 'get_1d_acceptance', 'get_horizontal_acceptance',
           'get_vertical_acceptance', 'get_momentum_acceptance']


def get_acceptance(ring, planes, npoints, amplitudes, nturns=1024,
                   refpts=None, dp=None, offset=None, bounds=None,
                   grid_mode=GridMode.RADIAL, use_mp=False, verbose=True,
                   start_method=None, divider=2, shift_zero=1.0e-9):
    """2D Acceptance calculation

    Computes the acceptance at ``repfts`` observation points.
    Scalar inputs can be used for 1D grid. The grid can be
    changed using grid_mode input:

    - ``at.GridMode.CARTESIAN``: :math:`(x,y)` grid
    - ``at.GridMode.RADIAL``: :math:`(r,\\theta)` grid
    - ``at.GridMode.RECURSIVE``: :math:`(r,\\theta)` recursive boundary search

    Example:
            >>> bf,sf,gf = ring.get_acceptance(planes, npoints, amplitudes)
            >>> plt.plot(*gf,'.')
            >>> plt.plot(*sf,'.')
            >>> plt.plot(*bf)
            >>> plt.show()

    args:
        ring:            ring use for tracking
        planes:          max. dimension 2, defines the plane where to search
                         for the acceptance allowed values are
                         ``x,xp,y,yp,dp,ct``
        npoints:         number of points in each dimension shape
                         ``(len(planes),)``
        amplitudes:      max. amplitude for ``RADIAL`` and ``CARTESIAN`` or
                         initial step
                         in ``RECURSIVE`` in each dimension. Shape is
                         ``(len(planes),)`` for ``RADIAL/RECURSIVE`` grid.
                         ``amplitude = sqrt(x**2+y**2)``

    keyword args:
        nturns=1024:     Number of turns for the tracking
        refpts=None:     Observation refpts, default start of the machine
        dp=None:         static momentum offset
        offset=None:     initial orbit, default closed orbit
        bounds=None:     defines the tracked range ``range=bounds*amplitude``,
                         it can be use to select quadrants for example
                         default values are
                         ``GridMode.CARTESIAN=((-1,1),(0,1))``
                         ``GridMode.RADIAL/RECURSIVE=((0,1),(pi,0))``
        grid_mode:       ``at.GridMode.CARTESIAN/RADIAL`` track full vector
                         (default). ``at.GridMode.RECURSIVE``: recursive search
        use_mp=False:    Use python multiprocessing (``patpass``, default use
                         ``lattice_pass``). In case multi-processing is not
                         enabled ``GridMode`` is forced to
                         ``RECURSIVE`` (most efficient in single core)
        divider=2:       Value of the divider used in ``RECURSIVE`` boundary
                         search
        verbose=True:    Print out some inform
        start_method:    This parameter allows to change the python
                         multiprocessing start method, ``default=None`` uses
                         the python defaults that is considered safe.
                         Available parameters are ``fork, spawn, forkserver``.
                         Default for linux is ``fork``, default for MacOS and
                         Windows is ``spawn``. ``fork`` may used for MacOS to
                         speed-up the calculation or to solve Runtime Errors,
                         however it is considered unsafe.

    return:
        The units depend on the selected planes and are the same as for the 6D
        particle coordinates. 3 elements tuples containing:

        - **boundary**: numpy array with shape ``(2,n)`` containing the 2D
          acceptance
        - **survived**: numpy array with shape ``(2,n)`` containing the initial
          coordinates of the particles that survived
        - **grid**: numpy array with shape ``(2,n)`` containing the initial
          coordinate of all the particles that

        In case ``len(refpts)>1`` a list of array of length ``len(refpts)``
        is returned.
    """
    kwargs = {}
    if not use_mp:
        grid_mode = GridMode.RECURSIVE
    elif start_method is not None:
        kwargs['start_method'] = start_method

    if verbose:
        nproc = multiprocessing.cpu_count()
        print('\n{0} cpu found for acceptance calculation'.format(nproc))
        if use_mp:
            nprocu = nproc
            print('Multi-process acceptance calculation selected...')
            if nproc == 1:
                print('Consider use_mp=False for single core computations')
        else:
            nprocu = 1
            print('Single process acceptance calculation selected...')
            if nproc > 1:
                print('Consider use_mp=True for parallelized computations')
        np = numpy.atleast_1d(npoints)
        na = 2
        if len(np) == 2:
            na = np[1]
        npp = numpy.prod(npoints)
        rpp = 2*numpy.ceil(numpy.log2(np[0]))*numpy.ceil(na/nprocu)
        mpp = npp/nprocu
        if rpp > mpp:
            cond = (grid_mode is GridMode.RADIAL or
                    grid_mode is GridMode.CARTESIAN)
        else:
            cond = grid_mode is GridMode.RECURSIVE
        if rpp > mpp and not cond:
            print('The estimated load for grid mode is {0}'.format(mpp))
            print('The estimated load for recursive mode is {0}'.format(rpp))
            print('{0} or {1} is recommended'.format(GridMode.RADIAL,
                                                     GridMode.CARTESIAN))
        elif rpp < mpp and not cond:
            print('The estimated load for grid mode is {0}'.format(mpp))
            print('The estimated load for recursive mode is {0}'.format(rpp))
            print('{0} is recommended'.format(GridMode.RECURSIVE))

    boundary = []
    survived = []
    grid = []
    if refpts is not None:
        rp = ring.uint32_refpts(refpts)
    else:
        rp = numpy.atleast_1d(refpts)
    for r in rp:
        b, s, g = boundary_search(ring, planes, npoints, amplitudes,
                                  nturns=nturns, obspt=r, dp=dp,
                                  offset=offset, bounds=bounds,
                                  grid_mode=grid_mode, use_mp=use_mp,
                                  verbose=verbose, divider=divider,
                                  shift_zero=shift_zero, **kwargs)
        boundary.append(b)
        survived.append(s)
        grid.append(g)
    if len(rp) == 1:
        return boundary[0], survived[0], grid[0]
    else:
        return boundary, survived, grid


def get_1d_acceptance(ring, plane, resolution, amplitude, nturns=1024, dp=None,
                      refpts=None, grid_mode=GridMode.RADIAL, use_mp=False,
                      verbose=False, start_method=None, divider=2):
    """1D Acceptance calculation

    Computes the acceptance at ``repfts`` observation points.
    Scalar inputs required.

    See also :meth:`at.acceptance.acceptance.get_acceptance`

    args:
        ring:            ring use for tracking
        plane:           defines the axis where to search for the acceptance
                         allowed values are ``x,xp,y,yp,dp,ct``
        resolution:      minimum distance between 2 grid points
        amplitude:       max. amplitude for ``RADIAL`` and ``CARTESIAN`` or
                         initial step in ``RECURSIVE``

    keyword args:
        nturns=1024:     Number of turns for the tracking
        refpts=None:     Observation refpts, default start of the machine
        dp=None:         static momentum offset
        grid_mode:       ``at.GridMode.CARTESIAN/RADIAL`` track full vector
                         (default). ``at.GridMode.RECURSIVE``: recursive search
        use_mp=False:    Use python multiprocessing (``patpass``, default use
                         ``lattice_pass``). In case multi-processing is not
                         enabled ``GridMode`` is forced to
                         ``RECURSIVE`` (most efficient in single core)
        divider=2:       Value of the divider used in ``RECURSIVE`` boundary
                         search
        verbose=True:    Print out some inform
        start_method:    This parameter allows to change the python
                         multiprocessing start method, ``default=None`` uses
                         the python defaults that is considered safe.
                         Available parameters are ``fork, spawn, forkserver``.
                         Default for linux is ``fork``, default for MacOS and
                         Windows is ``spawn``. ``fork`` may used for MacOS to
                         speed-up the calculation or to solve Runtime Errors,
                         however it is considered unsafe.

    return:
        The units depend on the selected planes and are the same as for the 6D
        particle coordinates. 3 elements tuples containing:

        - **boundary**: numpy array with shape ``(,n)`` containing the 1D
          acceptance
        - **survived**: numpy array with shape ``(,n)`` containing the initial
          coordinates of the particles that survived
        - **grid**: numpy array with shape ``(,n)`` containing the initial
          coordinate of all the particles that

        In case ``len(refpts)>1`` a list of array of length ``len(refpts)`` is
        returned. The boundary output is squeezed to an array with shape
        ``(len(refpts),2)``.
    """

    assert len(numpy.atleast_1d(plane)) == 1, \
        '1D acceptance: single plane required'
    assert numpy.isscalar(resolution), '1D acceptance: scalar args required'
    assert numpy.isscalar(amplitude), '1D acceptance: scalar args required'
    npoint = numpy.ceil(amplitude/resolution)
    b, s, g = get_acceptance(ring, plane, npoint, amplitude,
                             nturns=nturns, dp=dp, refpts=refpts,
                             grid_mode=grid_mode, use_mp=use_mp,
                             verbose=verbose, start_method=start_method,
                             divider=2, shift_zero=0.0)
    return numpy.squeeze(b), s, g


def get_horizontal_acceptance(ring, *args, **kwargs):
    """Horizontal acceptance calculation

    Computes the acceptance at ``repfts`` observation points.
    Scalar inputs required.

    See also :meth:`at.acceptance.acceptance.get_acceptance`
    and :meth:`at.acceptance.acceptance.get_1d_acceptance`

    args:
        ring:            ring use for tracking
        resolution:      minimum distance between 2 grid points
        amplitude:       max. amplitude for ``RADIAL`` and ``CARTESIAN`` or
                         initial step in ``RECURSIVE``

    keyword args:
        nturns=1024:     Number of turns for the tracking
        refpts=None:     Observation refpts, default start of the machine
        dp=None:         static momentum offset
        grid_mode:       ``at.GridMode.CARTESIAN/RADIAL`` track full vector
                         (default). ``at.GridMode.RECURSIVE``: recursive search
        use_mp=False:    Use python multiprocessing (``patpass``, default use
                         ``lattice_pass``). In case multi-processing is not
                         enabled ``GridMode`` is forced to
                         ``RECURSIVE`` (most efficient in single core)
        divider=2:       Value of the divider used in ``RECURSIVE`` boundary
                         search
        verbose=True:    Print out some inform
        start_method:    This parameter allows to change the python
                         multiprocessing start method, ``default=None`` uses
                         the python defaults that is considered safe.
                         Available parameters are ``fork, spawn, forkserver``.
                         Default for linux is ``fork``, default for MacOS and
                         Windows is ``spawn``. ``fork`` may used for MacOS to
                         speed-up the calculation or to solve Runtime Errors,
                         however it is considered unsafe.

    return:
        3 elements tuples containing:

        - **boundary**: numpy array with shape ``(,n)`` containing the 1D
          acceptance
        - **survived**: numpy array with shape ``(,n)`` containing the initial
          coordinates of the particles that survived
        - **grid**: numpy array with shape ``(,n)`` containing the initial
          coordinate of all the particles that

        In case ``len(refpts)>1`` a list of array of length ``len(refpts)`` is
        returned. The boundary output is squeezed to an array with shape
        ``(len(refpts),2)``
    """
    return get_1d_acceptance(ring, 'x', *args, **kwargs)


def get_vertical_acceptance(ring, *args, **kwargs):
    """Vertical acceptance calculation

    Computes the acceptance at ``repfts`` observation points.
    Scalar inputs required.

    See also :meth:`at.acceptance.acceptance.get_acceptance`
    and :meth:`at.acceptance.acceptance.get_1d_acceptance`

    args:
        ring:            ring use for tracking
        resolution:      minimum distance between 2 grid points
        amplitude:       max. amplitude for ``RADIAL`` and ``CARTESIAN`` or
                         initial step in ``RECURSIVE``

    keyword args:
        nturns=1024:     Number of turns for the tracking
        refpts=None:     Observation refpts, default start of the machine
        dp=None:         static momentum offset
        grid_mode:       ``at.GridMode.CARTESIAN/RADIAL`` track full vector
                         (default). ``at.GridMode.RECURSIVE``: recursive search
        use_mp=False:    Use python multiprocessing (``patpass``, default use
                         ``lattice_pass``). In case multi-processing is not
                         enabled ``GridMode`` is forced to
                         ``RECURSIVE`` (most efficient in single core)
        divider=2:       Value of the divider used in ``RECURSIVE`` boundary
                         search
        verbose=True:    Print out some inform
        start_method:    This parameter allows to change the python
                         multiprocessing start method, ``default=None`` uses
                         the python defaults that is considered safe.
                         Available parameters are ``fork, spawn, forkserver``.
                         Default for linux is ``fork``, default for MacOS and
                         Windows is ``spawn``. ``fork`` may used for MacOS to
                         speed-up the calculation or to solve Runtime Errors,
                         however it is considered unsafe.

    return:
        3 elements tuples containing:

        - **boundary**: numpy array with shape ``(,n)`` containing the 1D
          acceptance
        - **survived**: numpy array with shape ``(,n)`` containing the initial
          coordinates of the particles that survived
        - **grid**: numpy array with shape ``(,n)`` containing the initial
          coordinate of all the particles that

        In case ``len(refpts)>1`` a list of array of length ``len(refpts)`` is
        returned. The boundary output is squeezed to an array with shape
        ``(len(refpts),2)``
    """
    return get_1d_acceptance(ring, 'y', *args, **kwargs)


def get_momentum_acceptance(ring, *args, **kwargs):
    """Momentum acceptance calculation

    Computes the acceptance at ``repfts`` observation points.
    Scalar inputs required.

    See also :meth:`at.acceptance.acceptance.get_acceptance`
    and :meth:`at.acceptance.acceptance.get_1d_acceptance`
    This method is used in :meth:`at.acceptance.touschek.get_lifetime`

    args:
        ring:            ring use for tracking
        resolution:      minimum distance between 2 grid points
        amplitude:       max. amplitude for ``RADIAL`` and ``CARTESIAN`` or
                         initial step in ``RECURSIVE``

    keyword args:
        nturns=1024:     Number of turns for the tracking
        refpts=None:     Observation refpts, default start of the machine
        dp=None:         static momentum offset
        grid_mode:       ``at.GridMode.CARTESIAN/RADIAL`` track full vector
                         (default). ``at.GridMode.RECURSIVE``: recursive search
        use_mp=False:    Use python multiprocessing (``patpass``, default use
                         ``lattice_pass``). In case multi-processing is not
                         enabled ``GridMode`` is forced to
                         ``RECURSIVE`` (most efficient in single core)
        divider=2:       Value of the divider used in ``RECURSIVE`` boundary
                         search
        verbose=True:    Print out some inform
        start_method:    This parameter allows to change the python
                         multiprocessing start method, ``default=None`` uses
                         the python defaults that is considered safe.
                         Available parameters are ``fork, spawn, forkserver``.
                         Default for linux is ``fork``, default for MacOS and
                         Windows is ``spawn``. ``fork`` may used for MacOS to
                         speed-up the calculation or to solve Runtime Errors,
                         however it is considered unsafe.

    return:
        3 elements tuples containing:

        - **boundary**: numpy array with shape ``(,n)`` containing the 1D
          acceptance
        - **survived**: numpy array with shape ``(,n)`` containing the initial
          coordinates of the particles that survived
        - **grid**: numpy array with shape ``(,n)`` containing the initial
          coordinate of all the particles that

        In case ``len(refpts)>1`` a list of array of length ``len(refpts)`` is
        returned. The boundary output is squeezed to an array with shape
        ``(len(refpts),2)``
    """
    return get_1d_acceptance(ring, 'dp', *args, **kwargs)


Lattice.get_acceptance = get_acceptance
Lattice.get_horizontal_acceptance = get_horizontal_acceptance
Lattice.get_vertical_acceptance = get_vertical_acceptance
Lattice.get_momentum_acceptance = get_momentum_acceptance
