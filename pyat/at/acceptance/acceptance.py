"""Acceptance computation"""
import numpy
from .boundary import GridMode
# noinspection PyProtectedMember
from .boundary import boundary_search
from typing import Optional, Sequence
import multiprocessing
from ..lattice import Lattice, Refpts, frequency_control, AtError


__all__ = ['get_acceptance', 'get_1d_acceptance', 'get_horizontal_acceptance',
           'get_vertical_acceptance', 'get_momentum_acceptance']


@frequency_control
def get_acceptance(
        ring: Lattice, planes, npoints, amplitudes,
        nturns: Optional[int] = 1024,
        refpts: Optional[Refpts] = None,
        dp: Optional[float] = None,
        offset: Sequence[float] = None, bounds=None,
        grid_mode: Optional[GridMode] = GridMode.RADIAL,
        use_mp: Optional[bool] = False,
        verbose: Optional[bool] = True,
        divider: Optional[int] = 2,
        shift_zero: Optional[float] = 1.0e-6,
        start_method: Optional[str] = None,
):
    # noinspection PyUnresolvedReferences
    r"""Computes the acceptance at ``repfts`` observation points

    Parameters:
        ring:           Lattice definition
        planes:         max. dimension 2, Plane(s) to scan for the acceptance.
          Allowed values are: ``'x'``, ``'xp'``, ``'y'``,
          ``'yp'``, ``'dp'``, ``'ct'``
        npoints:        (len(planes),) array: number of points in each
          dimension
        amplitudes:     (len(planes),) array: set the search range:

          * :py:attr:`GridMode.CARTESIAN/RADIAL <.GridMode.RADIAL>`:
            max. amplitude
          * :py:attr:`.GridMode.RECURSIVE`: initial step
        nturns:         Number of turns for the tracking
        refpts:         Observation points. Default: start of the machine
        dp:             static momentum offset
        offset:         initial orbit. Default: closed orbit
        bounds:         defines the tracked range: range=bounds*amplitude.
          It can be use to select quadrants. For example, default values are:

          * :py:attr:`.GridMode.CARTESIAN`: ((-1, 1), (0, 1))
          * :py:attr:`GridMode.RADIAL/RECURSIVE <.GridMode.RADIAL>`: ((0, 1),
            (:math:`\pi`, 0))
        grid_mode:      defines the evaluation grid:

          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp:         Use python multiprocessing (:py:func:`.patpass`,
          default use :py:func:`.lattice_pass`).
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for MacOS and
          Windows is ``'spawn'``. ``'fork'`` may used for MacOS to speed-up
          the calculation or to solve runtime errors, however  it is
          considered unsafe.

    Returns:
        boundary:   (2,n) array: 2D acceptance
        survived:   (2,n) array: Coordinates of surviving particles
        tracked:    (2,n) array: Coordinates of tracked particles

    In case of multiple refpts, return values are lists of arrays, with one
    array per ref. point.

    Examples:

        >>> bf,sf,gf = ring.get_acceptance(planes, npoints, amplitudes)
        >>> plt.plot(*gf,'.')
        >>> plt.plot(*sf,'.')
        >>> plt.plot(*bf)
        >>> plt.show()

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
    """
    kwargs = {}
    if start_method is not None:
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
    if offset is not None:
        try:
            offset = numpy.broadcast_to(offset, (len(rp), 6))
        except ValueError:
            msg = ('offset and refpts have incoherent '
                   'shapes: {0}, {1}'.format(numpy.shape(offset),
                                             numpy.shape(refpts)))
            raise AtError(msg)
    else:
        offset=[None for _ in rp]
    for r, o in zip(rp, offset):
        b, s, g = boundary_search(ring, planes, npoints, amplitudes,
                                  nturns=nturns, obspt=r, dp=dp,
                                  offset=o, bounds=bounds,
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


def get_1d_acceptance(
        ring: Lattice, plane: str, resolution: float, amplitude: float,
        nturns: Optional[int] = 1024,
        refpts: Optional[Refpts] = None,
        dp: Optional[float] = None,
        offset: Sequence[float] = None,
        grid_mode: Optional[GridMode] = GridMode.RADIAL,
        use_mp: Optional[bool] = False,
        verbose: Optional[bool] = False,
        divider: Optional[int] = 2,
        shift_zero: Optional[float] = 1.0e-6,
        start_method: Optional[str] = None,

):
    r"""Computes the 1D acceptance at ``refpts`` observation points

    See :py:func:`get_acceptance`

    Parameters:
        ring:           Lattice definition
        plane:          Plane to scan for the acceptance.
          Allowed values are: ``'x'``, ``'xp'``, ``'y'``, ``'yp'``, ``'dp'``,
          ``'ct'``
        resolution:     Minimum distance between 2 grid points
        amplitude:      Search range:

          * :py:attr:`GridMode.CARTESIAN/RADIAL <.GridMode.RADIAL>`:
            max. amplitude
          * :py:attr:`.GridMode.RECURSIVE`: initial step
        nturns:         Number of turns for the tracking
        refpts:         Observation points. Default: start of the machine
        dp:             static momentum offset
        offset:         initial orbit. Default: closed orbit
        grid_mode:      defines the evaluation grid:

          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp:         Use python multiprocessing (:py:func:`.patpass`,
          default use :py:func:`.lattice_pass`). In case multi-processing
          is not enabled, ``grid_mode`` is forced to
          :py:attr:`.GridMode.RECURSIVE` (most efficient in single core)
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for MacOS and
          Windows is ``'spawn'``. ``'fork'`` may used for MacOS to speed-up
          the calculation or to solve runtime errors, however  it is considered
          unsafe.

    Returns:
        boundary:   (len(refpts),2) array: 1D acceptance
        tracked:    (n,) array: Coordinates of tracked particles
        survived:   (n,) array: Coordinates of surviving particles

    In case of multiple ``tracked`` and ``survived`` are lists of arrays,
    with one array per ref. point.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
    """
    if not use_mp:
        grid_mode = GridMode.RECURSIVE
    assert len(numpy.atleast_1d(plane)) == 1, \
        '1D acceptance: single plane required'
    assert numpy.isscalar(resolution), '1D acceptance: scalar args required'
    assert numpy.isscalar(amplitude), '1D acceptance: scalar args required'
    npoint = numpy.ceil(amplitude/resolution)
    if grid_mode is not GridMode.RECURSIVE:
        assert npoint > 1, \
            'Grid has only one point: increase amplitude or reduce resolution'
    b, s, g = get_acceptance(ring, plane, npoint, amplitude,
                             nturns=nturns, dp=dp, refpts=refpts,
                             grid_mode=grid_mode, use_mp=use_mp,
                             verbose=verbose, start_method=start_method,
                             divider=divider, shift_zero=shift_zero,
                             offset=offset)
    return numpy.squeeze(b), s, g


def get_horizontal_acceptance(ring: Lattice,
                              resolution: float, amplitude: float,
                              *args, **kwargs):
    r"""Computes the 1D horizontal acceptance at ``refpts`` observation points

    See :py:func:`get_acceptance`

    Parameters:
        ring:           Lattice definition
        resolution:     Minimum distance between 2 grid points
        amplitude:      Search range:

          * :py:attr:`GridMode.CARTESIAN/RADIAL <.GridMode.RADIAL>`:
            max. amplitude
          * :py:attr:`.GridMode.RECURSIVE`: initial step

    Keyword Args:
        nturns:         Number of turns for the tracking
        refpts:         Observation points. Default: start of the machine
        dp:             static momentum offset
        offset:         initial orbit. Default: closed orbit
        grid_mode:      defines the evaluation grid:

          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp:         Use python multiprocessing (:py:func:`.patpass`,
          default use :py:func:`.lattice_pass`). In case multi-processing
          is not enabled, ``grid_mode`` is forced to
          :py:attr:`.GridMode.RECURSIVE` (most efficient in single core)
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for MacOS and
          Windows is ``'spawn'``. ``'fork'`` may used for MacOS to speed-up
          the calculation or to solve runtime errors, however  it is considered
          unsafe.

    Returns:
        boundary:   (len(refpts),2) array: 1D acceptance
        tracked:    (n,) array: Coordinates of tracked particles
        survived:   (n,) array: Coordinates of surviving particles

    In case of multiple ``tracked`` and ``survived`` are lists of arrays,
    with one array per ref. point.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
    """
    return get_1d_acceptance(ring, 'x', resolution, amplitude, *args, **kwargs)


def get_vertical_acceptance(ring: Lattice,
                            resolution: float, amplitude: float,
                            *args, **kwargs):
    r"""Computes the 1D vertical acceptance at refpts observation points

    See :py:func:`get_acceptance`

    Parameters:
        ring:           Lattice definition
        resolution:     Minimum distance between 2 grid points
        amplitude:      Search range:

          * :py:attr:`GridMode.CARTESIAN/RADIAL <.GridMode.RADIAL>`:
            max. amplitude
          * :py:attr:`.GridMode.RECURSIVE`: initial step

    Keyword Args:
        nturns:         Number of turns for the tracking
        refpts:         Observation points. Default: start of the machine
        dp:             static momentum offset
        offset:         initial orbit. Default: closed orbit
        grid_mode:      defines the evaluation grid:

          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp:         Use python multiprocessing (:py:func:`.patpass`,
          default use :py:func:`.lattice_pass`). In case multi-processing
          is not enabled, ``grid_mode`` is forced to
          :py:attr:`.GridMode.RECURSIVE` (most efficient in single core)
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for MacOS and
          Windows is ``'spawn'``. ``'fork'`` may used for MacOS to speed-up
          the calculation or to solve runtime errors, however  it is considered
          unsafe.

    Returns:
        boundary:   (len(refpts),2) array: 1D acceptance
        tracked:    (n,) array: Coordinates of tracked particles
        survived:   (n,) array: Coordinates of surviving particles

    In case of multiple ``tracked`` and ``survived`` are lists of arrays,
    with one array per ref. point.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
    """
    return get_1d_acceptance(ring, 'y', resolution, amplitude, *args, **kwargs)


def get_momentum_acceptance(ring: Lattice,
                            resolution: float, amplitude: float,
                            *args, **kwargs):
    r"""Computes the 1D momentum acceptance at refpts observation points

    See :py:func:`get_acceptance`

    Parameters:
        ring:           Lattice definition
        resolution:     Minimum distance between 2 grid points
        amplitude:      Search range:

          * :py:attr:`GridMode.CARTESIAN/RADIAL <.GridMode.RADIAL>`:
            max. amplitude
          * :py:attr:`.GridMode.RECURSIVE`: initial step

    Keyword Args:
        nturns:         Number of turns for the tracking
        refpts:         Observation points. Default: start of the machine
        dp:             static momentum offset
        offset:         initial orbit. Default: closed orbit
        grid_mode:      defines the evaluation grid:

          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp:         Use python multiprocessing (:py:func:`.patpass`,
          default use :py:func:`.lattice_pass`). In case multi-processing is
          not enabled, ``grid_mode`` is forced to
          :py:attr:`.GridMode.RECURSIVE` (most efficient in single core)
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for MacOS and
          Windows is ``'spawn'``. ``'fork'`` may used for MacOS to speed-up
          the calculation or to solve runtime errors, however  it is considered
          unsafe.

    Returns:
        boundary:   (len(refpts),2) array: 1D acceptance
        tracked:    (n,) array: Coordinates of tracked particles
        survived:   (n,) array: Coordinates of surviving particles

    In case of multiple ``tracked`` and ``survived`` are lists of arrays,
    with one array per ref. point.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
    """
    return get_1d_acceptance(ring, 'dp', resolution, amplitude,
                             *args, **kwargs)


Lattice.get_acceptance = get_acceptance
Lattice.get_horizontal_acceptance = get_horizontal_acceptance
Lattice.get_vertical_acceptance = get_vertical_acceptance
Lattice.get_momentum_acceptance = get_momentum_acceptance
