"""Acceptance computation"""

from __future__ import annotations

__all__ = [
    "get_acceptance",
    "get_1d_acceptance",
    "get_horizontal_acceptance",
    "get_vertical_acceptance",
    "get_momentum_acceptance",
]

import multiprocessing
from typing import Sequence

import numpy as np

from .boundary import GridMode
from ..tracking import MPMode, gpu_core_count

# noinspection PyProtectedMember
from .boundary import boundary_search
from ..lattice import Lattice, Refpts, frequency_control


@frequency_control
def get_acceptance(
    ring: Lattice,
    planes,
    npoints,
    amplitudes,
    nturns: int = 1024,
    refpts: Refpts | None = None,
    dp: float | None = None,
    offset: Sequence[float] | None = None,
    bounds=None,
    grid_mode: GridMode = GridMode.RADIAL,
    use_mp: bool | MPMode = False,
    gpu_pool: list[int] | None = None,
    verbose: bool = False,
    divider: int = 2,
    shift_zero: float = 1.0e-6,
    start_method: str | None = None,
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
          It can be used to select quadrants. For example, default values are:

          * :py:attr:`.GridMode.CARTESIAN`: ((-1, 1), (0, 1))
          * :py:attr:`GridMode.RADIAL/RECURSIVE <.GridMode.RADIAL>`: ((0, 1),
            (:math:`\pi`, 0))
        grid_mode:      defines the evaluation grid:

          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp:         Flag to activate CPU or GPU multiprocessing
          (default: False)
        gpu_pool:       List of GPU id to use when use_mp is
          :py:attr:`at.tracking.MPMode.GPU`. If None specified, if gets
          first GPU.
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for macOS and
          Windows is ``'spawn'``. ``'fork'`` may be used on macOS to speed up
          the calculation or to solve runtime errors, however it is
          considered unsafe.

    Returns:
        boundary:   (2,n) array: 2D acceptance
        survived:   (2,n) array: Coordinates of surviving particles
        tracked:    (2,n) array: Coordinates of tracked particles

    In case of multiple refpts, return values are lists of arrays, with one
    array per ref. point.

    Examples:

        >>> bf, sf, gf = ring.get_acceptance(planes, npoints, amplitudes)
        >>> plt.plot(*gf, ".")
        >>> plt.plot(*sf, ".")
        >>> plt.plot(*bf)
        >>> plt.show()

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
       * When multiple ``refpts`` are provided particles are first
         projected to the beginning of the ring with tracking. Then,
         all particles are tracked up to ``nturns``. This allows to
         do most of the work in a single function call and allows for
         full parallelization.
    """
    kwargs = {}
    if start_method is not None:
        kwargs["start_method"] = start_method

    # For backward compatibility (use_mp can be a boolean)
    if use_mp is True:
        use_mp = MPMode.CPU

    if verbose:
        nproc = multiprocessing.cpu_count()
        print(f"\n{nproc} cpu found for acceptance calculation")
        if use_mp is MPMode.CPU:
            nprocu = nproc
            print("Multi-process acceptance calculation selected...")
            if nproc == 1:
                print("Consider use_mp=False for single core computations")
        elif use_mp is MPMode.GPU:
            nprocu = gpu_core_count(gpu_pool)
            print(f"\n{nprocu} GPU cores found")
            print("GPU acceptance calculation selected...")
            kwargs["gpu_pool"] = gpu_pool if gpu_pool is not None else [0]
        else:
            nprocu = 1
            print("Single process acceptance calculation selected...")
            if nproc > 1:
                print("Consider use_mp=True for parallelized computations")
        npts = np.atleast_1d(npoints)
        na = 2
        if len(npts) == 2:
            na = npts[1]
        npp = np.prod(npoints)
        rpp = 2 * np.ceil(np.log2(npts[0])) * np.ceil(na / nprocu)
        mpp = npp / nprocu
        if rpp > mpp:
            cond = grid_mode is GridMode.RADIAL or grid_mode is GridMode.CARTESIAN
        else:
            cond = grid_mode is GridMode.RECURSIVE
        if rpp > mpp and not cond:
            print(f"The estimated load for grid mode is {mpp}")
            print(f"The estimated load for recursive mode is {rpp}")
            print(f"{GridMode.RADIAL} or {GridMode.CARTESIAN} is recommended")
        elif rpp < mpp and not cond:
            print(f"The estimated load for grid mode is {mpp}")
            print(f"The estimated load for recursive mode is {rpp}")
            print(f"{GridMode.RECURSIVE} is recommended")

    b, s, g = boundary_search(
        ring,
        planes,
        npoints,
        amplitudes,
        nturns=nturns,
        obspt=refpts,
        dp=dp,
        offset=offset,
        bounds=bounds,
        grid_mode=grid_mode,
        use_mp=use_mp,
        verbose=verbose,
        divider=divider,
        shift_zero=shift_zero,
        **kwargs,
    )
    return b, s, g


def get_1d_acceptance(
    ring: Lattice,
    plane: str,
    resolution: float,
    amplitude: float,
    nturns: int | None = 1024,
    refpts: Refpts | None = None,
    dp: float | None = None,
    offset: Sequence[float] = None,
    grid_mode: GridMode | None = GridMode.RADIAL,
    use_mp: bool | MPMode = False,
    gpu_pool: list[int] | None = None,
    verbose: bool | None = False,
    divider: int | None = 2,
    shift_zero: float | None = 1.0e-6,
    start_method: str | None = None,
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
        use_mp:         Flag to activate CPU or GPU multiprocessing
          (default: False). In case multiprocessing is not enabled,
          ``grid_mode`` is forced to :py:attr:`.GridMode.RECURSIVE`
          (most efficient in single core)
        gpu_pool:       List of GPU id to use when use_mp is
          :py:attr:`at.tracking.MPMode.GPU`. If None specified, if gets
          first GPU.
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for macOS and
          Windows is ``'spawn'``. ``'fork'`` may be used on macOS to speed up
          the calculation or to solve runtime errors, however it is considered
          unsafe.

    Returns:
        boundary:   (len(refpts),2) array: 1D acceptance
        survived:   (n,) array: Coordinates of surviving particles
        tracked:    (n,) array: Coordinates of tracked particles

    In case of multiple ``tracked`` and ``survived`` are lists of arrays,
    with one array per ref. point.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
       * When multiple ``refpts`` are provided particles are first
         projected to the beginning of the ring with tracking. Then,
         all particles are tracked up to ``nturns``. This allows to
         do most of the work in a single function call and allows for
         full parallelization.
    """
    if not use_mp:
        if verbose:
            print("No parallel calculation selected, force to GridMode.RECURSIVE")
        grid_mode = GridMode.RECURSIVE
    assert len(np.atleast_1d(plane)) == 1, "1D acceptance: single plane required"
    assert np.isscalar(resolution), "1D acceptance: scalar args required"
    assert np.isscalar(amplitude), "1D acceptance: scalar args required"
    npoint = np.ceil(amplitude / resolution)
    if grid_mode is not GridMode.RECURSIVE:
        assert npoint > 1, (
            "Grid has only one point: increase amplitude or reduce resolution"
        )
    b, s, g = get_acceptance(
        ring,
        plane,
        npoint,
        amplitude,
        nturns=nturns,
        dp=dp,
        refpts=refpts,
        grid_mode=grid_mode,
        use_mp=use_mp,
        gpu_pool=gpu_pool,
        verbose=verbose,
        start_method=start_method,
        divider=divider,
        shift_zero=shift_zero,
        offset=offset,
    )
    return np.squeeze(b), s, g


def get_horizontal_acceptance(
    ring: Lattice, resolution: float, amplitude: float, *args, **kwargs
):
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
          default use :py:func:`.lattice_pass`). In case multiprocessing
          is not enabled, ``grid_mode`` is forced to
          :py:attr:`.GridMode.RECURSIVE` (most efficient in single core)
        gpu_pool:       List of GPU id to use when use_mp is
          :py:attr:`at.tracking.MPMode.GPU`. If None specified, if gets
          first GPU.
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for macOS and
          Windows is ``'spawn'``. ``'fork'`` may be used on macOS to speed up
          the calculation or to solve runtime errors, however it is considered
          unsafe.

    Returns:
        boundary:   (len(refpts),2) array: 1D acceptance
        survived:   (n,) array: Coordinates of surviving particles
        tracked:    (n,) array: Coordinates of tracked particles

    In case of multiple ``tracked`` and ``survived`` are lists of arrays,
    with one array per ref. point.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
       * When multiple ``refpts`` are provided particles are first
         projected to the beginning of the ring with tracking. Then,
         all particles are tracked up to ``nturns``. This allows to
         do most of the work in a single function call and allows for
         full parallelization.
    """
    return get_1d_acceptance(ring, "x", resolution, amplitude, *args, **kwargs)


def get_vertical_acceptance(
    ring: Lattice, resolution: float, amplitude: float, *args, **kwargs
):
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
          default use :py:func:`.lattice_pass`). In case multiprocessing
          is not enabled, ``grid_mode`` is forced to
          :py:attr:`.GridMode.RECURSIVE` (most efficient in single core)
        gpu_pool:       List of GPU id to use when use_mp is
          :py:attr:`at.tracking.MPMode.GPU`. If None specified, if gets
          first GPU.
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for macOS and
          Windows is ``'spawn'``. ``'fork'`` may be used on macOS to speed up
          the calculation or to solve runtime errors, however it is considered
          unsafe.

    Returns:
        boundary:   (len(refpts),2) array: 1D acceptance
        survived:   (n,) array: Coordinates of surviving particles
        tracked:    (n,) array: Coordinates of tracked particles

    In case of multiple ``tracked`` and ``survived`` are lists of arrays,
    with one array per ref. point.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
       * When multiple ``refpts`` are provided particles are first
         projected to the beginning of the ring with tracking. Then,
         all particles are tracked up to ``nturns``. This allows to
         do most of the work in a single function call and allows for
         full parallelization.
    """
    return get_1d_acceptance(ring, "y", resolution, amplitude, *args, **kwargs)


def get_momentum_acceptance(
    ring: Lattice, resolution: float, amplitude: float, *args, **kwargs
):
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
          default use :py:func:`.lattice_pass`). In case multiprocessing is
          not enabled, ``grid_mode`` is forced to
          :py:attr:`.GridMode.RECURSIVE` (most efficient in single core)
        gpu_pool:       List of GPU id to use when use_mp is
          :py:attr:`at.tracking.MPMode.GPU`. If None specified, if gets
          first GPU.
        verbose:        Print out some information
        divider:        Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for macOS and
          Windows is ``'spawn'``. ``'fork'`` may be used on macOS to speed up
          the calculation or to solve runtime errors, however it is considered
          unsafe.

    Returns:
        boundary:   (len(refpts),2) array: 1D acceptance
        survived:   (n,) array: Coordinates of surviving particles
        tracked:    (n,) array: Coordinates of tracked particles

    In case of multiple ``tracked`` and ``survived`` are lists of arrays,
    with one array per ref. point.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
       * When multiple ``refpts`` are provided particles are first
         projected to the beginning of the ring with tracking. Then,
         all particles are tracked up to ``nturns``. This allows to
         do most of the work in a single function call and allows for
         full parallelization.
    """
    return get_1d_acceptance(ring, "dp", resolution, amplitude, *args, **kwargs)


Lattice.get_acceptance = get_acceptance
Lattice.get_horizontal_acceptance = get_horizontal_acceptance
Lattice.get_vertical_acceptance = get_vertical_acceptance
Lattice.get_momentum_acceptance = get_momentum_acceptance
