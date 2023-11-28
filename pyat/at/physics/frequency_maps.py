"""
Frequency analysis (FMAP) using .harmonic_analysis lib
and pyat parallel tracking (patpass)
"""

# orblancog
# generates the frequency and diffusion map for a given ring
# 2023jan16 tracking is parallel (patpass), frequency analysis is serial
# 2022jun07 serial version

from at.tracking import patpass
from at.physics import find_orbit
import numpy
from warnings import warn
from at.lattice import AtWarning


# Jaime Coello de Portugal (JCdP) frequency analysis implementation
from .harmonic_analysis import get_tunes_harmonic

__all__ = ['fmap_parallel_track']


def fmap_parallel_track(ring,
                        coords=[-10, 10, -10, 10],
                        steps=[100, 100],
                        scale='linear',
                        turns=512,
                        orbit=None,
                        add_offset6D=numpy.zeros(6),
                        verbose=False,
                        lossmap=False,
                        **kwargs
                        ):
    r"""Computes frequency maps

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
    Additionally, a numpy array (*add_offset6D*) with 6 values could be used to
    arbitrarily offset the initial coordinates of every particle.

    A dictionary with particle losses is saved for every vertical offset.
    See the :py:func:`.patpass` documentation.

    Parameters:
        ring:     a valid pyat ring
        coords:   default [-10,10,-10,10] in mm
        steps:    (xsteps, ysteps): number of steps in each plane
        scale:    default 'linear'; or 'non-linear'
        turns:    default 512
        orbit:    If :py:obj:`None`, the closed orbit is computed and added to
          the coordinates
        add_offset6D: default numpy.zeros((6,1))
        verbose:  prints additional info
        lossmap:  default false
    Optional:
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
    kwargs.setdefault('pool_size', None)

    # https://github.com/atcollab/at/pull/608
    if 'ncpu' in kwargs:
        warn(AtWarning('ncpu argument is deprecated; use pool_size instead'))
        kwargs['pool_size'] = kwargs.pop('ncpu')

    if orbit is None:
        # closed orbit values are not returned. It seems not necessary here
        orbit, _ = find_orbit(ring)
        print(f'Closed orbit:\t{orbit}')

    # verboseprint to check flag only once
    verboseprint = print if verbose else lambda *a, **k: None

    # tns is the variable used in the frequency analysis
    # turns is the input variable from user
    # nturns is twice tns in order to get the tune in the first and second
    #     part of the tracking
    tns = turns
    nturns = 2*tns

    # scale the inumpyut coordinates
    xscale = 1e-3
    yscale = 1e-3

    # returned array
    xy_nuxy_lognudiff_array = numpy.empty([])
    loss_map_array = numpy.empty([])

    # verify steps
    if numpy.count_nonzero(steps) != 2:
        raise ValueError('steps can not be zero')

    xmin = numpy.minimum(coords[0], coords[1])
    xmax = numpy.maximum(coords[0], coords[1])
    ymin = numpy.minimum(coords[2], coords[3])
    ymax = numpy.maximum(coords[2], coords[3])
    xsteps = steps[0]
    ysteps = steps[1]

    # define rectangle and x,y step size
    if scale == "nonlinear":
        verboseprint('non linear steps')
        xinterval = 1.0*(xmax - xmin)
        yinterval = 1.0*(ymax - ymin)
        xauxsteps = numpy.arange(-0.5, 0.0 + 1e-9, 1./xsteps)
        yauxsteps = numpy.arange(-0.5, 0.0 + 1e-9, 1./ysteps)
        xnrsteps = 10**xauxsteps
        xnrsteps = xnrsteps / sum(xnrsteps)
        ynrsteps = 10**yauxsteps
        ynrsteps = ynrsteps / sum(ynrsteps)
        xscalesteps = numpy.concatenate((xnrsteps,
                                        numpy.flip(xnrsteps)), axis=None)
        yscalesteps = numpy.concatenate((ynrsteps,
                                        numpy.flip(ynrsteps)), axis=None)
        ixarray = xmin + 0.5 * xinterval * numpy.cumsum(xscalesteps)
        iyarray = ymin + 0.5 * yinterval * numpy.cumsum(yscalesteps)
    else:  # scale == 'linear', ignore any other
        verboseprint('linear steps')
        # get the intervals
        xstep = 1.0*(xmax - xmin)/xsteps
        ystep = 1.0*(ymax - ymin)/ysteps
        ixarray = numpy.arange(xmin, xmax+1e-6, xstep)
        iyarray = numpy.arange(ymin, ymax+1e-6, ystep)
    lenixarray = len(ixarray)
    leniyarray = len(iyarray)

    print("Start tracking and frequency analysis")

    # tracking in parallel multiple x coordinates with the same y coordinate
    for iy, iy_index in zip(iyarray, range(leniyarray)):
        print(f'Tracked particles {abs(-100.0*iy_index/leniyarray):.1f} %')
        verboseprint("y =", iy)
        z0 = numpy.zeros((lenixarray, 6))  # transposed, and C-aligned
        z0 = z0 + add_offset6D + orbit
        # add 1 nm to tracking to avoid zeros in array for the ideal lattice
        z0[:, 0] = z0[:, 0] + xscale*ixarray + 1e-9
        z0[:, 2] = z0[:, 2] + yscale*iy + 1e-9

        verboseprint("tracking ...")
        # z0.T is Fortran-aligned
        if lossmap:
            # patpass output changes when losses flag is true
            zOUT, dictloss = \
                patpass(ring, z0.T, nturns, losses=lossmap, **kwargs)
            loss_map_array = numpy.append(loss_map_array, dictloss)
        else:
            zOUT = patpass(ring, z0.T, nturns, **kwargs)

        # start of serial frequency analysis
        for ix_index in range(lenixarray):  # cycle over the track results
            # check if nan in arrays
            array_sum = numpy.sum(zOUT[:, ix_index, 0])
            array_has_nan = numpy.isnan(array_sum)
            if array_has_nan:
                verboseprint("array has nan")
                continue

            # get one valid particle
            z1 = zOUT[:, ix_index, 0]

            # remove mean values
            # get the first turn in x
            xfirst = z1[0, 0:tns]
            xfirst = xfirst - numpy.mean(xfirst)
            # pxfirst = z1[1, 0:tns]
            # pxfirst = pxfirst - numpy.mean(pxfirst)
            # xfirstpart = xfirst + 1j*pxfirst
            # get the last turns in x
            xlast = z1[0, tns:2*tns]
            xlast = xlast - numpy.mean(xlast)
            # pxlast = z1[1, tns:2*tns]
            # pxlast = pxlast - numpy.mean(pxlast)
            # xlastpart = xlast + 1j*pxlast

            # get the first turn in y
            yfirst = z1[2, 0:tns]
            yfirst = yfirst - numpy.mean(yfirst)
            # pyfirst = z1[3, 0:tns]
            # pyfirst = pyfirst - numpy.mean(pyfirst)
            # yfirstpart = yfirst + 1j*pyfirst
            # get the last turns in y
            ylast = z1[2, tns:2*tns]
            ylast = ylast - numpy.mean(ylast)
            # pylast = z1[3, tns:2*tns]
            # pylast = pylast - numpy.mean(pylast)
            # ylastpart = ylast + 1j*pylast

            # calc frequency from array,
            # jump the cycle is no frequency is found
            xfreqfirst = get_tunes_harmonic(xfirst)
            if len(xfreqfirst) == 0:
                verboseprint("  No frequency")
                continue
            # xfreqlast  = PyNAFF.naff(xlastpart,tns,1,0,False)
            xfreqlast = get_tunes_harmonic(xlast)
            if len(xfreqlast) == 0:
                verboseprint("  No frequency")
                continue
            # yfreqfirst = PyNAFF.naff(yfirstpart,tns,1,0,False)
            yfreqfirst = get_tunes_harmonic(yfirst)
            if len(yfreqfirst) == 0:
                verboseprint("  No frequency")
                continue
            # yfreqlast  = PyNAFF.naff(ylastpart,tns,1,0,False)
            yfreqlast = get_tunes_harmonic(ylast)
            if len(yfreqlast) == 0:
                verboseprint("  No frequency")
                continue
            verboseprint("NAFF results, (x,y)=",
                         ixarray[ix_index],
                         iy,
                         "\nH freq. first part =\t", xfreqfirst[0],
                         "\nH freq. last part =\t", xfreqlast[0],
                         "\nV freq. first part =\t", yfreqfirst[0],
                         "\nV freq. last part =\t", yfreqlast[0])

            # metric
            xdiff = xfreqlast[0] - xfreqfirst[0]
            ydiff = yfreqlast[0] - yfreqfirst[0]
            nudiff = 0.5*numpy.log10((xdiff*xdiff + ydiff*ydiff)/tns)
            # min max diff
            if nudiff > -2:
                nudiff = -2
            if nudiff < -10:
                nudiff = -10
            # save diff
            xy_nuxy_lognudiff_array = numpy.append(xy_nuxy_lognudiff_array,
                                                   [ixarray[ix_index],
                                                    iy,
                                                    xfreqfirst[0],
                                                    yfreqfirst[0],
                                                    xdiff,
                                                    ydiff,
                                                    nudiff])

    # first element is garbage
    xy_nuxy_lognudiff_array = numpy.delete(xy_nuxy_lognudiff_array, 0)
    # reshape for plots and output files
    xy_nuxy_lognudiff_array = xy_nuxy_lognudiff_array.reshape(-1, 7)

    return xy_nuxy_lognudiff_array, loss_map_array
# the end
