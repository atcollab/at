""" Frequency analysis (FMAP). """


# orblancog
# generates the frequency and diffusion map for a given ring
# 2024apr16 create get_freqmap
# 2023jan16 tracking is parallel (patpass), frequency analysis is serial
# 2022jun07 serial version
# 2024jun21 adds get_freqmap


from __future__ import annotations

import multiprocessing
import time
from typing import Optional, Sequence
from warnings import warn

import numpy

from ..acceptance.boundary import (GridMode, get_parts, get_plane_index,
                                   get_survived, grid_configuration,
                                   set_ring_orbit)
from ..lattice import AtError, AtWarning, Lattice, Refpts
from ..physics import find_orbit
from ..tracking import patpass
# Jaime Coello de Portugal (JCdP) frequency analysis implementation
from .harmonic_analysis import get_tunes_harmonic

_pdict = {"x": 0, "xp": 1, "y": 2, "yp": 3, "dp": 4, "ct": 5}


__all__ = ["fmap_parallel_track", "get_freqmap"]


def get_freqmap(
    ring: Lattice,
    planes: list,
    npoints: numpy.ndarray,
    amplitudes: numpy.ndarray,
    bounds: any = None,
    nturns: Optional[int] = 512,
    deltap: Optional[float] = None,
    offset: Sequence[float] = None,
    refpts: Optional[Refpts] = 0,
    grid_mode: Optional[GridMode] = GridMode.CARTESIAN,
    use_mp: Optional[bool] = True,
    verbose: Optional[bool] = False,
    shift_zero: Optional[float] = 0.0e-6,
    start_method: Optional[str] = None,
):
    # noinspection PyUnresolvedReferences
    r"""Computes the frequency map at ``repfts`` observation points.

    Parameters:
        ring:           Lattice definition
        planes:         max. dimension 2, Plane(s) to scan for the acceptance.
          Allowed values are: ``'x'``, ``'xp'``, ``'y'``,
          ``'yp'``, ``'dp'``, ``'ct'``
          e.g. ['x','y']
        npoints:        (len(planes),) array: number of points in each
          dimension
          e.g. numpy.array([200,200])
        amplitudes:     (len(planes),) array: set the search range:
          * :py:attr:`GridMode.CARTESIAN/RADIAL <.GridMode.RADIAL>`:
            max. amplitude
          * :py:attr:`.GridMode.RECURSIVE`: initial step
          e.g. numpy.array([10e-3, 10e-3])
    Keyword arguments:
        nturns:         Number of turns for the tracking. Default 512
        refpts:         Observation points. Default: start of the machine
        deltap:             static momentum offset
        offset:         (len(refpts),6) array: initial orbit. Default: closed orbit
        bounds:         defines the tracked range: range=bounds*amplitude.
          It can be use to select quadrants. For example, default values are:
          * :py:attr:`.GridMode.CARTESIAN`: ((-1, 1), (0, 1))
          * :py:attr:`GridMode.RADIAL/RECURSIVE <.GridMode.RADIAL>`: ((0, 1),
            (:math:`\pi`, 0))
        grid_mode:      defines the evaluation grid:
          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp:         Default True.
          Use python multiprocessing (:py:func:`.patpass`,
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
        fmaout:  (len(refpts)) list: containing the result of every reference point.
            Per every reference point there is a tuple containing,
              1) (npoints*npoints, 7) array with the frequency map,
                 each row has
                  `[xcoor, ycoor, nux, ny, dnux, dnuy,
                   log10(sqrt(sum(dnu**2)/turns)) ]`
              2) the grid
              3) the dictionary of particle losses


    Examples:

        >>> fmaout = at.get_freqmap(ring,
                     ['x','y'],
                     numpy.array([200,200]),
                     numpy.array([10e-3,10e-3]),
                     nturns=512,
                     refpts = numpy.array([0,1,2]),
                     bounds=[[-1,1],[-1,1]],
                     verbose=True,
                     )
        >>> fmadata_at_start = fmaout[0][0]

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
    """
    kwargs = {}
    if start_method is not None:
        kwargs["start_method"] = start_method

    if verbose:
        nproc = multiprocessing.cpu_count()
        print(f"\n{0} cpu found for acceptance calculation".format(nproc))
        if use_mp:
            print("Multi-process acceptance calculation selected...")
            if nproc == 1:
                print("Consider use_mp=False for single core computations")
        else:
            print("Single process acceptance calculation selected...")
            if nproc > 1:
                print("Consider use_mp=True for parallelized computations")

    survived = []
    grid = []
    if refpts is not None:
        rps = ring.uint32_refpts(refpts)
    else:
        rps = numpy.atleast_1d(refpts)
    if offset is not None:
        try:
            offset = numpy.broadcast_to(offset, (len(rps), 6))
        except ValueError as incoherent_shape:
            msg = f"offset and refpts have incoherent shapes: {0}, {1}".format(
                numpy.shape(offset), numpy.shape(refpts)
            )
            raise AtError(msg) from incoherent_shape
    else:
        _, offset = find_orbit(ring, refpts)
        planesi = numpy.atleast_1d(get_plane_index(planes))
        offset[:, planesi[0]] = offset[:, planesi[0]] + 1e-9
        offset[:, planesi[1]] = offset[:, planesi[1]] + 1e-9
    dataobs = []
    t00 = time.time()
    for obspt, offset0 in zip(rps, offset):
        offset, newring = set_ring_orbit(ring, deltap, obspt, offset0)
        config = grid_configuration(
            planes, npoints, amplitudes, grid_mode, bounds=bounds, shift_zero=shift_zero
        )
        if verbose:
            print("\nRunning grid frequency search:")
            if obspt is None:
                print(f"Element {0}, obspt={1}".format(ring[0].FamName, 0))
            else:
                print(f"Element {0}, obspt={1}".format(ring[obspt].FamName, obspt))
            print(f"The grid mode is {0}".format(config.mode))
            print(f"The planes are {0}".format(config.planes))
            print(f"Number of steps are {0}".format(config.shape))
            print(f"The maximum amplitudes are {0}".format(config.amplitudes))
            print(f"The maximum boundaries are {0}".format(config.bounds))
            print(f"Number of turns is {0}".format(nturns))
            print(f"The initial offset is {0} with deltap={1}".format(offset, deltap))
        parts, grid = get_parts(config, offset)
        rout, _, tdl = ring.track(
            parts, nturns=2 * nturns, losses=True, use_mp=use_mp, **kwargs
        )

        mask = get_survived(parts, newring, nturns, use_mp, **kwargs)
        survived = grid.grid[:, mask]

        planesi = numpy.atleast_1d(get_plane_index(planes))
        tunes = numpy.zeros((len(planesi), 2, len(numpy.where(mask)[0])))

        for i, theplane in enumerate(planesi):
            tunes[i, 0] = get_tunes_harmonic(
                rout[theplane, mask, :, 0:nturns], use_mp=use_mp, **kwargs
            )
            tunes[i, 1] = get_tunes_harmonic(
                rout[theplane, mask, :, nturns : 2 * nturns], use_mp=use_mp, **kwargs
            )
        # metric
        diffplane1 = tunes[0, 0, :] - tunes[0, 1, :]
        diffplane2 = tunes[1, 0, :] - tunes[1, 1, :]
        nudiff = 0.5 * numpy.log10(
            (diffplane1 * diffplane1 + diffplane2 * diffplane2) / nturns
        )
        # set min-max
        nudiff = numpy.clip(nudiff, -10, -2)

        # now return rout,tp,tdl,grid,tunes
        firstturns = 0
        dataobs.append(
            (
                numpy.vstack(
                    (
                        survived,
                        tunes[0, firstturns],
                        tunes[1, firstturns],
                        diffplane1,
                        diffplane2,
                        nudiff,
                    )
                ).T,
                grid,
                tdl,
            )
        )
    if verbose:
        print(f"Calculation took {0}".format(time.time() - t00))

    return dataobs


def fmap_parallel_track(
    ring: Lattice,
    coords=[-10, 10, -10, 10],
    steps: list = [100, 100],
    scale: str = "linear",
    turns: int = 512,
    orbit: numpy.ndarray = None,
    add_offset6D: numpy.ndarray = numpy.zeros(6),
    add_offset6d: numpy.ndarray = numpy.zeros(6),
    verbose: bool = False,
    lossmap: bool = False,
    **kwargs: dict[str, any],
) -> numpy.ndarray:
    r"""Computes frequency maps.
    This function calculates the norm of the transverse tune variation per turn
    for a particle tracked along a ring with a set of offsets in the initial
    coordinates.

    It returns a numpy array containing 7 columns:

    | ``[xoffset, yoffset, nux, nuy, dnux, dnuy,``
    | ``log10(sqrt(dnux*dnux + dnuy*dnuy)/tns )]``

    for every tracked particle that survives 2*turns.
    Particles lost before 2*turns are ignored.
    The log10 tune variation is limited to the interval from -10 to -2;
    values above or below are set to the closer limit.

    The transverse offsets are given inside a rectangular coordinate window

    ``coords=[xmin, xmax, ymin, *ymax]``

    in millimeters, where the window is divided in n steps=[xsteps,ysteps].
    For each yoffset, particle tracking over the whole set of xoffsets is done
    in parallel.

    The frequency analysis uses a library that is not parallelized.

    The closed orbit is calculated and added to the
    initial particle offset of every particle. Otherwise, one could set
    ``orbit = numpy.zeros(6)`` to avoid it.
    Additionally, a numpy array (*add_offset6d*) with 6 values could be used to
    arbitrarily offset the initial coordinates of every particle.

    A dictionary with particle losses is saved for every vertical offset.
    See the :py:func:`.patpass` documentation.

    Parameters:
        ring:     a valid pyat ring
    Optional:
        coords:   default [-10,10,-10,10] in mm
        steps:    (xsteps, ysteps): number of steps in each plane
        scale:    default 'linear'; or 'non-linear'
        turns:    default 512
        orbit:    If :py:obj:`None`, the closed orbit is computed and added to
          the coordinates
        add_offset6d: default numpy.zeros((6,1))
        verbose:  prints additional info
        lossmap:  default false
        pool_size:number of processes. See :py:func:`.patpass`

    Returns:
        xy_nuxy_lognudiff_array: numpy array with columns
         `[xcoor, ycoor, nux, ny, dnux, dnuy,
         log10(sqrt(sum(dnu**2)/turns)) ]`
        loss_map_array : experimental format.
          if *loss_map* is :py:obj:`True`, it returns the losses dictionary
          provided by :py:func:`.patpass` per every vertical offset.
          if *loss_map* is :py:obj:`False`, it returns a one-element list.

    .. warning::

       points with NaN tracking results or non-defined frequency in x,y
       are ignored.

    .. warning:: loss map format is experimental
    """

    # https://github.com/atcollab/at/pull/608
    kwargs.setdefault("pool_size", None)

    # https://github.com/atcollab/at/pull/608
    if "ncpu" in kwargs:
        warn(AtWarning("ncpu argument is deprecated; use pool_size instead"))
        kwargs["pool_size"] = kwargs.pop("ncpu")

    # deprecating add_offset6D because does not comply with the cammel name standard
    if "add_offset6D" in kwargs:
        warn(AtWarning("add_offset6D argument is deprecated; use add_offset6d instead"))
        kwargs["add_offset6d"] = kwargs.pop("add_offset6D")
    add_offset6d = kwargs.pop("add_offset6d", numpy.zeros(6))

    if orbit is None:
        # closed orbit values are not returned. It seems not necessary here
        orbit, _ = find_orbit(ring)
        print(f"Closed orbit:\t{orbit}")

    # verboseprint to check flag only once
    verboseprint = print if verbose else lambda *a, **k: None

    # tns is the variable used in the frequency analysis
    # turns is the input variable from user
    # nturns is twice tns in order to get the tune in the first and second
    #     part of the tracking
    tns = turns
    nturns = 2 * tns

    # scale the inumpyut coordinates
    xscale = 1e-3
    yscale = 1e-3

    # returned array
    xy_nuxy_lognudiff_array = numpy.empty([])
    loss_map_array = numpy.empty([])

    # verify steps
    if numpy.count_nonzero(steps) != 2:
        raise ValueError("steps can not be zero")

    xmin = numpy.minimum(coords[0], coords[1])
    xmax = numpy.maximum(coords[0], coords[1])
    ymin = numpy.minimum(coords[2], coords[3])
    ymax = numpy.maximum(coords[2], coords[3])
    xsteps = steps[0]
    ysteps = steps[1]

    # define rectangle and x,y step size
    if scale == "nonlinear":
        verboseprint("non linear steps")
        xinterval = 1.0 * (xmax - xmin)
        yinterval = 1.0 * (ymax - ymin)
        xauxsteps = numpy.arange(-0.5, 0.0 + 1e-9, 1.0 / xsteps)
        yauxsteps = numpy.arange(-0.5, 0.0 + 1e-9, 1.0 / ysteps)
        xnrsteps = 10**xauxsteps
        xnrsteps = xnrsteps / sum(xnrsteps)
        ynrsteps = 10**yauxsteps
        ynrsteps = ynrsteps / sum(ynrsteps)
        xscalesteps = numpy.concatenate((xnrsteps, numpy.flip(xnrsteps)), axis=None)
        yscalesteps = numpy.concatenate((ynrsteps, numpy.flip(ynrsteps)), axis=None)
        ixarray = xmin + 0.5 * xinterval * numpy.cumsum(xscalesteps)
        iyarray = ymin + 0.5 * yinterval * numpy.cumsum(yscalesteps)
    else:  # scale == 'linear', ignore any other
        verboseprint("linear steps")
        # get the intervals
        xstep = 1.0 * (xmax - xmin) / xsteps
        ystep = 1.0 * (ymax - ymin) / ysteps
        ixarray = numpy.arange(xmin, xmax + 1e-6, xstep)
        iyarray = numpy.arange(ymin, ymax + 1e-6, ystep)
    lenixarray = len(ixarray)
    leniyarray = len(iyarray)

    print("Start tracking and frequency analysis")

    # tracking in parallel multiple x coordinates with the same y coordinate
    for iiy, iy_index in zip(iyarray, range(leniyarray)):
        print(f"Tracked particles {abs(-100.0*iy_index/leniyarray):.1f} %")
        verboseprint("y =", iiy)
        z00 = numpy.zeros((lenixarray, 6))  # transposed, and C-aligned
        z00 = z00 + add_offset6d + orbit
        # add 1 nm to tracking to avoid zeros in array for the ideal lattice
        z00[:, 0] = z00[:, 0] + xscale * ixarray + 1e-9
        z00[:, 2] = z00[:, 2] + yscale * iiy + 1e-9

        verboseprint("tracking ...")
        # z00.T is Fortran-aligned
        if lossmap:
            # patpass output changes when losses flag is true
            zout, dictloss = patpass(ring, z00.T, nturns, losses=lossmap, **kwargs)
            loss_map_array = numpy.append(loss_map_array, dictloss)
        else:
            zout = patpass(ring, z00.T, nturns, **kwargs)

        # start of serial frequency analysis
        for ix_index in range(lenixarray):  # cycle over the track results
            # check if nan in arrays
            array_sum = numpy.sum(zout[:, ix_index, 0])
            array_has_nan = numpy.isnan(array_sum)
            if array_has_nan:
                verboseprint("array has nan")
                continue

            # get one valid particle
            z11 = zout[:, ix_index, 0]

            # remove mean values
            # get the first turn in x
            xfirst = z11[0, 0:tns]
            xfirst = xfirst - numpy.mean(xfirst)
            # get the last turns in x
            xlast = z11[0, tns : 2 * tns]
            xlast = xlast - numpy.mean(xlast)

            # get the first turn in y
            yfirst = z11[2, 0:tns]
            yfirst = yfirst - numpy.mean(yfirst)
            # get the last turns in y
            ylast = z11[2, tns : 2 * tns]
            ylast = ylast - numpy.mean(ylast)

            # calc frequency from array,
            # jump the cycle is no frequency is found
            xfreqfirst = get_tunes_harmonic(xfirst)
            if len(xfreqfirst) == 0:
                verboseprint("  No frequency")
                continue
            xfreqlast = get_tunes_harmonic(xlast)
            if len(xfreqlast) == 0:
                verboseprint("  No frequency")
                continue
            yfreqfirst = get_tunes_harmonic(yfirst)
            if len(yfreqfirst) == 0:
                verboseprint("  No frequency")
                continue
            yfreqlast = get_tunes_harmonic(ylast)
            if len(yfreqlast) == 0:
                verboseprint("  No frequency")
                continue
            verboseprint(
                "NAFF results, (x,y)=",
                ixarray[ix_index],
                iiy,
                "\nH freq. first part =\t",
                xfreqfirst[0],
                "\nH freq. last part =\t",
                xfreqlast[0],
                "\nV freq. first part =\t",
                yfreqfirst[0],
                "\nV freq. last part =\t",
                yfreqlast[0],
            )

            # metric
            xdiff = xfreqlast[0] - xfreqfirst[0]
            ydiff = yfreqlast[0] - yfreqfirst[0]
            nudiff = 0.5 * numpy.log10((xdiff * xdiff + ydiff * ydiff) / tns)
            # min max diff
            if nudiff > -2:
                nudiff = -2
            if nudiff < -10:
                nudiff = -10
            # save diff
            xy_nuxy_lognudiff_array = numpy.append(
                xy_nuxy_lognudiff_array,
                [
                    ixarray[ix_index],
                    iiy,
                    xfreqfirst[0],
                    yfreqfirst[0],
                    xdiff,
                    ydiff,
                    nudiff,
                ],
            )

    # first element is garbage
    xy_nuxy_lognudiff_array = numpy.delete(xy_nuxy_lognudiff_array, 0)
    # reshape for plots and output files
    xy_nuxy_lognudiff_array = xy_nuxy_lognudiff_array.reshape(-1, 7)

    return xy_nuxy_lognudiff_array, loss_map_array


# the end
