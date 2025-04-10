"""
Functions used to
calculate the loss boundary for different
grid definitions
"""

from at.lattice import Lattice, AtError, AtWarning, Refpts
from typing import Optional, Sequence
from enum import Enum
import numpy
from scipy.ndimage import binary_dilation, binary_opening
from collections import namedtuple
import time
import warnings

__all__ = ["GridMode"]

_pdict = {"x": 0, "xp": 1, "y": 2, "yp": 3, "dp": 4, "ct": 5}


class GridMode(Enum):
    """
    Grid definition for 2D acceptance boundary search
    """

    RADIAL = 0  #: full [:math:`\:r, \theta\:`] grid
    CARTESIAN = 1  #: full [:math:`\:x, y\:`] grid
    RECURSIVE = 2  #: radial recursive search


def grid_config(planes, amplitudes, npoints, bounds, grid_mode, shift_zero):
    """
    Returns an object that defines the grid configuration
    """
    bounds = numpy.atleast_2d(bounds)
    bounds = tuple(map(tuple, bounds))
    shape = numpy.array(npoints, dtype=numpy.int32)
    d = {
        "planes": numpy.atleast_1d(planes),
        "planesi": numpy.atleast_1d(get_plane_index(planes)),
        "amplitudes": numpy.atleast_1d(amplitudes),
        "shape": numpy.atleast_1d(shape),
        "bounds": bounds,
        "mode": grid_mode,
        "shift_zero": shift_zero,
    }
    return namedtuple("config", d.keys())(**d)


def grid(grid, offset):
    """
    Returns a grid object
    """
    d = {"grid": numpy.atleast_2d(grid), "offset": numpy.atleast_1d(offset)}
    return namedtuple("grid", d.keys())(**d)


def get_plane_index(planes):
    """
    Converts plane to particle coordinate index
    """
    planesi = numpy.array([], dtype=numpy.int32)
    for i, p in enumerate(numpy.atleast_1d(planes)):
        if isinstance(p, str):
            try:
                planesi = numpy.append(planesi, _pdict[p])
            except KeyError:
                raise AtError("Allowed values for plane are x,xp,y,yp,dp,ct")
        else:
            raise AtError("Allowed values for plane are x,xp,y,yp,dp,ct")
    return planesi


def set_ring_orbit(ring, dp, obspt, orbit):
    """
    Returns a ring starting at obspt and initial
    closed orbit
    """
    if obspt is not None:
        assert numpy.isscalar(obspt), "Scalar value needed for obspt"
        ring = ring.rotate(obspt)

    if orbit is None:
        orbit = ring.find_orbit(dp=dp)[0]

    return orbit, ring


def grid_configuration(
    planes, npoints, amplitudes, grid_mode, bounds=None, shift_zero=1.0e-6
):
    """
    Return a grid configuration based on user input parameters, the ordering
    is as follows: CARTESIAN: (x,y), RADIAL/RECURSIVE (r, theta).
    Scalar inputs can be used for 1D grid
    """
    ndims = len(numpy.atleast_1d(planes))
    if ndims > 2 or ndims == 0:
        raise AtError("planes can have 1 or 2 element (1D or 2D aperture)")
    elif ndims == 1 and grid_mode is GridMode.RADIAL:
        grid_mode = GridMode.CARTESIAN

    if numpy.shape(numpy.atleast_1d(npoints)) != (ndims,):
        raise AtError("npoints shape should be (len(planes),)")
    if numpy.shape(numpy.atleast_1d(amplitudes)) != (ndims,):
        raise AtError("amplitudes shape should be (len(planes),)")
    if numpy.shape(numpy.atleast_2d(bounds)) != (ndims, 2) and bounds is not None:
        raise AtError("bounds shape should be (len(planes),2)")

    if grid_mode is GridMode.RADIAL or grid_mode is GridMode.RECURSIVE:
        if bounds is None:
            bounds = numpy.array([[0, 1], [numpy.pi, 0]])
        bounds[0][bounds[0] == 0] = 1.0e-6
    elif grid_mode is GridMode.CARTESIAN:
        if bounds is None:
            bounds = numpy.array([[p - 1, 1] for p in range(ndims)])
    else:
        raise AtError("GridMode {0} undefined.".format(grid_mode))

    config = grid_config(planes, amplitudes, npoints, bounds, grid_mode, shift_zero)
    return config


def get_parts(config, offset):
    """
    Generate a (6,n) particles array based on the grid configuration
    and initial offset, returns a grid object containing the grid and
    the offset from which the array can be reconstructed
    """

    def get_part_grid_uniform(bnd, np, amp):
        x = [numpy.linspace(*b, n) * a for a, b, n in zip(amp, bnd, np)]
        try:
            g1, g2 = numpy.meshgrid(*x)
        except ValueError:
            g1, g2 = numpy.meshgrid(x, 0.0)
        return numpy.array([g1.flatten(), g2.flatten()])

    def get_part_grid_radial(bnd, np, amp):
        x = [numpy.linspace(*b, n) for b, n in zip(bnd, np)]
        g1, g2 = numpy.meshgrid(*x)
        g1r = amp[0] * numpy.cos(g2) * g1
        g2r = amp[1] * numpy.sin(g2) * g1
        return numpy.array([g1r.flatten(), g2r.flatten()])

    pind = config.planesi
    amp = config.amplitudes
    np = config.shape
    bnd = config.bounds
    gm = config.mode

    if gm is GridMode.CARTESIAN:
        g = get_part_grid_uniform(bnd, np, amp)
    elif gm is GridMode.RADIAL:
        g = get_part_grid_radial(bnd, np, amp)
    parts = numpy.zeros((6, numpy.prod(np)))
    parts[pind, :] = [g[i] for i in range(len(pind))]
    offset = numpy.array(offset) + config.shift_zero
    parts = (parts.T + offset).T
    return parts, grid(g, offset[pind])


def get_survived(parts, ring, nturns, use_mp, **kwargs):
    """
    Track a grid through the ring and extract survived particles
    """
    _, _, td = ring.track(
        parts,
        nturns=nturns,
        losses=True,
        use_mp=use_mp,
        refpts=None,
        in_place=True,
        **kwargs,
    )
    return numpy.invert(td["loss_map"].islost)


def get_grid_boundary(mask, grid, config):
    """
    Compute the boundary of survided particles array
    """

    def nearest_order(grid):
        #  keep only max r for each angle
        angle = numpy.arctan2(*grid)
        norm = numpy.linalg.norm(grid.T, axis=1)
        val, inv = numpy.unique(angle, return_inverse=True)
        gf = numpy.zeros((2, len(val)))
        for i, v in enumerate(val):
            inds = numpy.where(inv == i)[0]
            ni = norm[inds]
            nim = numpy.where(ni == numpy.amax(ni))[0]
            ind = inds[nim][0]
            gf[:, i] = grid[:, ind]
        #  first sort by angle
        idx = numpy.argsort(numpy.arctan2(*gf))
        gf = gf[:, idx]
        #  now sort by closest neighbour on normalized grid
        x, y = gf[0, :].copy(), gf[1, :].copy()
        dxmin = min(numpy.diff(numpy.unique(x)))
        dymin = min(numpy.diff(numpy.unique(y)))
        iorder = [0]
        for i in range(1, len(x)):
            xnow = x[iorder[-1]]
            ynow = y[iorder[-1]]
            dd = numpy.sqrt(((x - xnow) / dxmin) ** 2 + ((y - ynow) / dymin) ** 2)
            if i <= 3:
                ic = [j for j in numpy.argsort(dd) if j not in iorder]
            else:
                direction = numpy.sign(iorder[-1] - iorder[-2])
                ic = [
                    j
                    for j in numpy.argsort(dd)
                    if j not in iorder and numpy.sign(j - iorder[-1]) == direction
                ]
            if len(ic) > 0:
                iorder.append(ic[0])
        #  finally connect both ends if distance within unit square
        xnow = x[iorder[-1]]
        ynow = y[iorder[-1]]
        xs = x[iorder[0]]
        ys = y[iorder[0]]
        dd = numpy.sqrt(((xs - xnow) / dxmin) ** 2 + ((ys - ynow) / dymin) ** 2)
        if dd < 1.5:
            iorder.append(iorder[0])
        gf = gf[:, iorder]
        return gf

    def search_bnd(ma, sa):
        bnd = numpy.zeros((2, 1))
        for j, m in enumerate(ma):
            bnd = sa[:, j]
            if not m and j > 0:
                bnd = sa[:, j - 1]
                break
        return bnd

    def grid_boundary(mask, grid, config, nclean=2):
        cnt = numpy.flip(config.shape)
        if len(cnt) == 1:
            return vector_boundary(mask, grid)
        mask2d = numpy.reshape(mask.copy(), cnt)
        #  remove isolated points
        for i in range(nclean):
            bnd1 = binary_opening(mask2d, numpy.ones((2, 1)))
            bnd2 = binary_opening(mask2d, numpy.ones((1, 2)))
        bnd = numpy.logical_and(bnd1, bnd2)
        k = numpy.zeros((3, 3), dtype=int)
        k[1] = 1
        k[:, 1] = 1
        bnd = numpy.logical_and(binary_dilation(bnd == 0, border_value=1), bnd)
        bnd = grid.grid[:, bnd.reshape(mask.shape)]
        return nearest_order(bnd)

    def vector_boundary(mask, grid):
        g = grid.grid
        xp, xn = g[0] >= 0, g[0] <= 0
        bp = search_bnd(mask[xp], g[:, xp])
        bn = search_bnd(numpy.flip(mask[xn]), numpy.flip(g[:, xn], axis=1))
        return numpy.squeeze([bn[0], bp[0]])

    def radial_boundary(mask, grid):
        angles = numpy.round(numpy.arctan2(*grid.grid), decimals=8)
        angles, invi = numpy.unique(angles, return_inverse=True)
        bnd = numpy.zeros((2, len(angles)))
        for i in range(len(angles)):
            sa = grid.grid[:, invi == i]
            ma = mask[invi == i]
            bnd[:, i] = search_bnd(ma, sa)
        return bnd

    if not numpy.any(mask):
        msg = (
            "No particle survived, please check your grid "
            "or lattice. Acceptance set to [0.0, 0.0]."
        )
        warnings.warn(AtWarning(msg))
        cnt = numpy.flip(config.shape)
        if len(cnt) == 1:
            return numpy.zeros(2)
        else:
            return numpy.zeros((2, 1))

    if config.mode is GridMode.RADIAL:
        return radial_boundary(mask, grid)
    elif config.mode is GridMode.CARTESIAN:
        return grid_boundary(mask, grid, config)
    else:
        raise AtError("GridMode {0} undefined.".format(grid.mode))


def grid_boundary_search(
    ring,
    planes,
    npoints,
    amplitudes,
    nturns=1024,
    obspt=None,
    dp=None,
    offset=None,
    bounds=None,
    grid_mode=GridMode.RADIAL,
    use_mp=False,
    verbose=True,
    shift_zero=1.0e-9,
    **kwargs,
):
    """
    Search for the boundary by tracking a grid
    """
    config = grid_configuration(
        planes, npoints, amplitudes, grid_mode, bounds=bounds, shift_zero=shift_zero
    )

    if verbose:
        kwargs["verbose"] = verbose
        print("\nRunning grid boundary search:")
        if len(obspt) == 1:
            if obspt[0] is None:
                print("Element {0}, obspt={1}".format(ring[0].FamName, 0))
            else:
                print("Element {0}, obspt={1}".format(ring[obspt].FamName, obspt))
        else:
            print(
                "{0} Elements from {1}, obspt={2} to {3}, obspt={4}".format(
                    len(obspt),
                    ring[obspt[0]].FamName,
                    obspt[0],
                    ring[obspt[-1]].FamName,
                    obspt[-1],
                )
            )
        print("The grid mode is {0}".format(config.mode))
        print("The planes are {0}".format(config.planes))
        print("Number of steps are {0}".format(config.shape))
        print("The maximum amplitudes are {0}".format(config.amplitudes))
        print("The maximum boundaries are {0}".format(config.bounds))

    t0 = time.time()
    allparts = []
    grids = []
    offsets = []

    for i, obs, orbit in zip(numpy.arange(len(obspt)), obspt, offset):
        orbit, newring = set_ring_orbit(ring, dp, obs, orbit)
        parts, grid = get_parts(config, orbit)
        obs = 0 if obs is None else obs
        dpp = 0.0 if dp is None else dp
        if verbose:
            print(
                "\r{4}/{5}: Projecting obs=({0}, {1}) to the start of the ring, "
                "the initial offset is {2} with dp={3}".format(
                    ring[obs].FamName, obs, orbit, dpp, i + 1, len(obspt)
                )
            )
        newring[: len(ring) - obs].track(
            parts, use_mp=use_mp, in_place=True, refpts=None, **kwargs
        )
        allparts.append(parts)
        grids.append(grid)
        offsets.append(orbit)
    if verbose:
        print("Starting the multi-turn tracking...")
    allparts = numpy.concatenate(allparts, axis=1)
    mask = get_survived(allparts, ring, nturns, use_mp, **kwargs)
    mask = numpy.split(mask, len(grids))
    survived = [g.grid[:, m] for g, m in zip(grids, mask)]
    boundary = [get_grid_boundary(m, g, config) for g, m in zip(grids, mask)]
    grids = [g.grid for g in grids]
    if verbose:
        print("Calculation took {0}".format(time.time() - t0))
    if len(obspt) == 1:
        return boundary[0], survived[0], grids[0]
    else:
        return boundary, survived, grids


def recursive_boundary_search(
    ring,
    planes,
    npoints,
    amplitudes,
    nturns=1024,
    obspt=None,
    dp=None,
    offset=None,
    bounds=None,
    use_mp=False,
    divider=2,
    verbose=True,
    shift_zero=1.0e-9,
    **kwargs,
):
    """
    Recursively search for the boundary in a given plane and direction (angle)
    """

    def search_boundary(
        planesi, angles, rtol, rsteps, nturns, offset, use_mp, **kwargs
    ):

        ftol = min(rtol / rsteps)
        cs = numpy.squeeze([numpy.cos(angles), numpy.sin(angles)])
        cs = numpy.around(cs, decimals=9)
        fact = numpy.ones(len(angles))
        survived = numpy.full(len(angles), True)
        part = numpy.zeros((6, len(angles)))
        grid = numpy.array([])
        mask = numpy.array([])

        while numpy.any(survived):
            for i, pi in enumerate(planesi):
                part[pi, survived] += cs[i, survived] * rsteps[i] * fact[survived]
            istracked = numpy.array(
                [
                    not numpy.any([numpy.allclose(p, g, rtol=1.0e-9) for g in grid.T])
                    for p in part[planesi].T
                ]
            )
            survived = numpy.array(
                [
                    numpy.any([numpy.allclose(p, m, rtol=1.0e-9) for m in mask.T])
                    for p in part[planesi].T
                ]
            )
            pt = part[:, istracked]
            grid = numpy.hstack([grid, pt[planesi]]) if grid.size else pt[planesi]
            ptmp = (pt.T + offset).T
            survived[istracked] = get_survived(ptmp, newring, nturns, use_mp, **kwargs)
            pm = part[:, numpy.logical_and(istracked, survived)]
            mask = numpy.hstack([mask, pm[planesi]]) if mask.size else pm[planesi]
            for i in range(len(angles)):
                if not survived[i] and fact[i] > ftol:
                    deltas = cs[:, i] * rsteps[:] * min(1, 2 * fact[i])
                    if numpy.any(abs(deltas) > abs(part[planesi, i])):
                        part[planesi, i] = numpy.zeros(len(planesi))
                    else:
                        for j, pi in enumerate(planesi):
                            part[pi, i] -= deltas[j]
                    survived[i] = True
                    fact[i] *= 1 / divider

        for i, pi in enumerate(planesi):
            part[pi] -= cs[i] * rsteps[i] * fact

        p = numpy.squeeze(part[planesi])
        return p, mask, grid

    offset, newring = set_ring_orbit(ring, dp, obspt, offset)
    config = grid_configuration(
        planes,
        npoints,
        amplitudes,
        GridMode.RECURSIVE,
        bounds=bounds,
        shift_zero=shift_zero,
    )
    rtol = min(numpy.atleast_1d(config.amplitudes / config.shape))
    rstep = config.amplitudes
    if len(numpy.atleast_1d(config.shape)) == 2:
        angles = numpy.linspace(*config.bounds[1], config.shape[1])
    else:
        angles = numpy.linspace(*config.bounds[1], 2)
    angles = numpy.atleast_1d(angles)

    if verbose:
        print("\nRunning recursive boundary search:")
        if obspt is None:
            print("Element {0}, obspt={1}".format(ring[0].FamName, 0))
        else:
            print("Element {0}, obspt={1}".format(ring[obspt].FamName, obspt))
        print("The grid mode is {0}".format(config.mode))
        print("The planes are {0}".format(config.planes))
        print(
            "Number of angles is {0} from {1} to {2} rad".format(
                len(angles), angles[0], angles[-1]
            )
        )
        print("The resolution of the search is {0}".format(rtol))
        print("The initial step size is {0}".format(rstep))
        print("The initial offset is {0} with dp={1}".format(offset, dp))

    t0 = time.time()
    result = search_boundary(
        config.planesi, angles, rtol, rstep, nturns, offset, use_mp, **kwargs
    )
    if verbose:
        print("Calculation took {0}".format(time.time() - t0))
    return result


def boundary_search(
    ring: Lattice,
    planes,
    npoints,
    amplitudes,
    nturns: Optional[int] = 1024,
    obspt: Optional[Refpts] = None,
    dp: Optional[float] = None,
    offset: Sequence[float] = None,
    bounds=None,
    grid_mode: Optional[GridMode] = GridMode.RADIAL,
    use_mp: Optional[bool] = False,
    verbose: Optional[bool] = True,
    shift_zero: Optional[float] = 1.0e-9,
    **kwargs,
):
    """
    Computes the loss boundary at a single point in the machine
    """
    if obspt is not None:
        rp = ring.uint32_refpts(obspt)
    else:
        rp = numpy.atleast_1d(obspt)
    if offset is not None:
        try:
            offset = numpy.broadcast_to(offset, (len(rp), 6))
        except ValueError:
            msg = "offset and refpts have incoherent " "shapes: {0}, {1}".format(
                numpy.shape(offset), numpy.shape(obspt)
            )
            raise AtError(msg)
    else:
        offset = [None for _ in rp]

    divider = kwargs.pop("divider", 2)
    if grid_mode is GridMode.RECURSIVE:
        boundary = []
        survived = []
        grid = []
        for r, o in zip(rp, offset):
            b, s, g = recursive_boundary_search(
                ring,
                planes,
                npoints,
                amplitudes,
                nturns=nturns,
                obspt=r,
                dp=dp,
                offset=o,
                bounds=bounds,
                use_mp=use_mp,
                verbose=verbose,
                divider=divider,
                shift_zero=shift_zero,
                **kwargs,
            )
            boundary.append(b)
            survived.append(s)
            grid.append(g)
        if len(rp) == 1:
            result = (boundary[0], survived[0], grid[0])
        else:
            result = (boundary, survived, grid)
    else:
        result = grid_boundary_search(
            ring,
            planes,
            npoints,
            amplitudes,
            nturns=nturns,
            obspt=rp,
            dp=dp,
            offset=offset,
            bounds=bounds,
            grid_mode=grid_mode,
            use_mp=use_mp,
            verbose=verbose,
            shift_zero=shift_zero,
            **kwargs,
        )
    return result
