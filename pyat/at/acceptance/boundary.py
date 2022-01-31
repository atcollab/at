import at
from at.lattice import AtError, AtWarning
from enum import Enum
import numpy
import warnings
from scipy.ndimage.morphology import binary_dilation, binary_opening
from collections import namedtuple
import time

__all__ = ['GridMode']

_pdict = {'x': 0, 'xp': 1,
          'y': 2, 'yp': 3,
          'dp': 4, 'ct': 5}


class GridMode(Enum):
    """"
    Class to defined the grid mode use when searching
    for the boundary
    """
    RADIAL = 0
    CARTESIAN = 1
    RECURSIVE = 2


def grid_config(planes, amplitudes, npoints, bounds, grid_mode):
    """"
    Returns an object that defines the grid configuration
    """
    bounds = numpy.atleast_2d(bounds)
    bounds = tuple(map(tuple, bounds))
    shape = numpy.array(npoints, dtype=numpy.int32)
    d = {'planes': numpy.atleast_1d(planes),
         'planesi': numpy.atleast_1d(get_plane_index(planes)),
         'amplitudes': numpy.atleast_1d(amplitudes),
         'shape': numpy.atleast_1d(shape),
         'bounds': bounds,
         'mode': grid_mode}
    return namedtuple('config', d.keys())(**d)


def grid(grid, offset):
    """"
    Returns a grid object
    """
    d = {'grid': numpy.atleast_2d(grid),
         'offset': numpy.atleast_1d(offset)}
    return namedtuple('grid', d.keys())(**d)


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
                raise AtError('Allowed values for plane are x,xp,y,yp,dp,ct')
        else:
            raise AtError('Allowed values for plane are x,xp,y,yp,dp,ct')
    return planesi


def set_ring_orbit(ring, dp, refpts, orbit):
    """
    Returns a ring starting at refpts and initial
    closed orbit
    """
    newring = ring.set_rf_frequency(dp=dp, copy=True)
    if refpts is not None:
        assert numpy.isscalar(refpts), 'Scalar value needed for refpts'
        newring = newring.rotate(refpts)
    if orbit is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            orbit, _ = ring.find_orbit(dp=dp)
    return orbit, newring


def grid_configuration(planes, npoints, amplitudes, grid_mode, bounds=None):
    """
    Return a grid configuration based on user input parameters, the ordering
    is as follows: CARTESIAN: (x,y), RADIAL/RECURSIVE (r, theta). Scalar inputs can
    be used for 1D grid
    """
    ndims = len(numpy.atleast_1d(planes))

    if numpy.shape(numpy.atleast_1d(npoints)) != (ndims,):
        raise AtError('npoints shape should be (len(planes),)')
    if numpy.shape(numpy.atleast_1d(amplitudes)) != (ndims,):
        raise AtError('amplitudes shape should be (len(planes),)')
    if (numpy.shape(numpy.atleast_2d(bounds)) != (ndims, 2)
       and bounds is not None):
        raise AtError('bounds shape should be (len(planes),2)')

    if ndims > 2 or ndims == 0:
        raise AtError('planes can have 1 or 2 element (1D or 2D aperture)')
    elif ndims == 1 and grid_mode is GridMode.RADIAL:
        grid_mode = GridMode.CARTESIAN

    if grid_mode is GridMode.RADIAL or grid_mode is GridMode.RECURSIVE:
        if bounds is None:
            bounds = numpy.array([[0, 1], [numpy.pi, 0]])
        bounds[0][bounds[0] == 0] = 1.0e-6
    elif grid_mode is GridMode.CARTESIAN:
        if bounds is None:
            bounds = numpy.array([[p-1, 1] for p in range(ndims)])
    else:
        raise AtError('GridMode {0} undefined.'.format(grid_mode))

    config = grid_config(planes, amplitudes, npoints,
                         bounds, grid_mode)
    return config


def get_parts(config, offset):
    """
    Generate a (6,n) particles array based on the grid configuration
    and initial offset, returns a grid object containing the grid and
    the offset from which the array can be reconstructed
    """
    def get_part_grid_uniform(bnd, np, amp):
        x = [numpy.linspace(*b, n)*a for a, b, n in zip(amp, bnd, np)]
        try:
            g1, g2 = numpy.meshgrid(*x)
        except ValueError:
            g1, g2 = numpy.meshgrid(x, 0.0)
        return numpy.array([g1.flatten(), g2.flatten()])

    def get_part_grid_radial(bnd, np, amp):
        x = [numpy.linspace(*b, n) for b, n in zip(bnd, np)]
        g1, g2 = numpy.meshgrid(*x)
        g1r = amp[0]*numpy.cos(g2)*g1
        g2r = amp[1]*numpy.sin(g2)*g1
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
    parts = (parts.T+offset).T
    return parts, grid(g, offset[pind])


def get_survived(parts, ring, nturns, method):
    """
    Track a grid through the ring and extract survived particles
    """
    if not (method is at.patpass or method is at.lattice_pass):
        raise AtError('Only patpass (multi-process) and lattice_pass '
                      '(single process) allowed for tracking method')
    pout = numpy.squeeze(method(ring, parts, nturns=nturns))
    if pout.ndim == 2:
        return numpy.invert(numpy.isnan(pout[0, -1]))
    else:
        return numpy.invert(numpy.isnan(pout[0, :, -1]))


def get_grid_boundary(mask, grid, config):
    """
    Compute the boundary of survided particles array
    """
    def nearest_order(grid):
        #  first sort by angle
        idx = numpy.argsort(numpy.arctan2(*grid))
        grid = grid[:, idx]
        #  now sort by closest neighbour on normalized grid
        x, y = grid[0, :].copy(), grid[1, :].copy()
        dxmin = min(numpy.diff(numpy.unique(x)))
        dymin = min(numpy.diff(numpy.unique(y)))
        iorder = [0]
        for i in range(1, len(x)):
            xnow = x[iorder[-1]]
            ynow = y[iorder[-1]]
            dd = numpy.sqrt(((x-xnow)/dxmin)**2+((y-ynow)/dymin)**2)
            ic = [j for j in numpy.argsort(dd) if j not in iorder]
            #  take only points inside unit square
            if dd[ic[0]] < numpy.sqrt(2)*1.1:
                iorder.append(ic[0])
        return grid[:, iorder]

    def search_bnd(ma, sa):
        bnd = numpy.zeros((2, 1))
        for j, m in enumerate(ma):
            bnd = sa[:, j]
            if not m and j > 0:
                bnd = sa[:, j-1]
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
        bnd = numpy.logical_and(binary_dilation(bnd == 0), bnd)
        bnd = grid.grid[:, bnd.reshape(mask.shape)]
        return nearest_order(bnd)

    def vector_boundary(mask, grid):
        g = grid.grid
        xp, xn = g[0] >= 0, g[0] <= 0
        bp = search_bnd(mask[xp], g[:, xp])
        bn = search_bnd(numpy.flip(mask[xn]),
                        numpy.flip(g[:, xn], axis=1))
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

    if config.mode is GridMode.RADIAL:
        return radial_boundary(mask, grid)
    elif config.mode is GridMode.CARTESIAN:
        return grid_boundary(mask, grid, config)
    else:
        raise AtError('GridMode {0} undefined.'.format(grid.mode))


def grid_boundary_search(ring, planes, npoints, amplitudes, nturns=1024,
                         refpts=None, dp=0, offset=None, bounds=None,
                         grid_mode=GridMode.RADIAL, use_mp=False,
                         verbose=True):
    """
    Search for the boundary by tracking a grid
    """
    if use_mp:
        method = at.patpass
    else:
        method = at.lattice_pass

    offset, newring = set_ring_orbit(ring, dp, refpts,
                                     offset)
    config = grid_configuration(planes, npoints, amplitudes,
                                grid_mode, bounds=bounds)

    if verbose:
        print('\nRunning grid boundary search:')
        if refpts is None:
            print('Element {0}, refpts={1}'.format(ring[0].FamName, 0))
        else:
            print('Element {0}, refpts={1}'.format(ring[refpts].FamName,
                                                   refpts))
        print('The grid mode is {0}'.format(config.mode))
        print('The planes are {0}'.format(config.planes))
        print('Number of steps are {0}'.format(config.shape))
        print('The maximum amplitudes are {0}'.format(config.amplitudes))
        print('The maximum boundaries are {0}'.format(config.bounds))
        print('The initial offset is {0} with dp={1}'.format(offset, dp))

    t0 = time.time()
    parts, grid = get_parts(config, offset)
    mask = get_survived(parts, newring, nturns, method)
    survived = grid.grid[:, mask]
    boundary = get_grid_boundary(mask, grid, config)
    if verbose:
        print('Calculation took {0}'.format(time.time()-t0))
    return boundary, survived, grid.grid


def recursive_boundary_search(ring, planes, npoints, amplitudes, nturns=1024,
                              refpts=None, dp=0, offset=None, bounds=None,
                              use_mp=False, verbose=True):
    """
    Recursively search for the boundary in a given plane and direction (angle)
    """
    def get_r(arr):
        r = numpy.linalg.norm(arr)
        return numpy.around(r, decimals=9)

    def search_boundary(ring, planesi, angles, rtol, rstep, nturns,
                        offset, method):

        rstep[rstep < rtol] = rtol
        if len(planesi) == 1:
            cs = numpy.atleast_2d(numpy.cos(angles)).T
            rsteps = numpy.array([numpy.squeeze(rstep) for c in cs])
        else:
            cs = numpy.squeeze([numpy.cos(angles), numpy.sin(angles)]).T
            rsteps = numpy.array([numpy.sqrt((rstep[0]*c[0])**2 +
                                             (rstep[1]*c[1])**2) for c in cs])
        cs = numpy.around(cs, decimals=9)

        steps = [{} for a in angles]
        survived = numpy.full(len(angles), True)
        istracked = numpy.full(len(angles), True)
        part = numpy.zeros((6, len(angles)))
        while numpy.any(survived):
            for pi, c in zip(planesi, cs.T):
                part[pi, survived] += c[survived]*rsteps[survived]
            for i in range(len(angles)):
                try:
                    survived[i] = steps[i][get_r(part[planesi, i])]
                    istracked[i] = False
                except KeyError:
                    istracked[i] = True

            ptmp = (part[:, istracked].T + offset).T
            survived[istracked] = get_survived(ptmp, newring, nturns, method)

            for i in range(len(angles)):
                rp = get_r(part[planesi, i])
                steps[i][rp] = survived[i]
                if not survived[i] and rsteps[i] > rtol:
                    part[planesi, i] -= cs[i]*min(2.0*rsteps[i], rp)
                    rsteps[i] *= 0.5
                    survived[i] = True

        for pi, c in zip(planesi, cs.T):
            part[pi] -= c*rsteps

        grid = []
        mask = []

        for a, stp in zip(angles, steps):
            for k, v in stp.items():
                x = numpy.around(k*numpy.cos(a), decimals=9)
                y = numpy.around(k*numpy.sin(a), decimals=9)
                if v:
                    mask.append([x, y])
                grid.append([x, y])
        p = numpy.squeeze(part[planesi])
        m = numpy.array(mask).T
        g = numpy.array(grid).T
        return p, m, g

    if use_mp:
        method = at.patpass
    else:
        method = at.lattice_pass
    offset, newring = set_ring_orbit(ring, dp, refpts, offset)
    config = grid_configuration(planes, npoints, amplitudes,
                                GridMode.RECURSIVE, bounds=bounds)
    rtol = min(numpy.atleast_1d(config.amplitudes/config.shape))
    rstep = config.amplitudes
    if len(numpy.atleast_1d(config.shape)) == 2:
        angles = numpy.linspace(*config.bounds[1], config.shape[1])
    else:
        angles = numpy.linspace(*config.bounds[1], 2)
    angles = numpy.atleast_1d(angles)

    if verbose:
        print('\nRunning recursive boundary search:')
        if refpts is None:
            print('Element {0}, refpts={1}'.format(ring[0].FamName, 0))
        else:
            print('Element {0}, refpts={1}'.format(ring[refpts].FamName,
                                                   refpts))
        print('The grid mode is {0}'.format(config.mode))
        print('The planes are {0}'.format(config.planes))
        print('Number of angles is {0} from {1} to {2} rad'.format(len(angles),
              angles[0], angles[-1]))
        print('The resolution of the search is {0}'.format(rtol))
        print('The initial step size is {0}'.format(rstep))
        print('The initial offset is {0} with dp={1}'.format(offset, dp))

    t0 = time.time()
    result = search_boundary(ring, config.planesi, angles, rtol, rstep,
                             nturns, offset, method)
    if verbose:
        print('Calculation took {0}'.format(time.time()-t0))
    return result


def boundary_search(ring, planes, npoints, amplitudes, nturns=1024,
                    refpts=None, dp=0, offset=None, bounds=None,
                    grid_mode=GridMode.RADIAL, use_mp=False, verbose=True):
    """
    Computes the loss boundary at a single point in the machine
    """
    if grid_mode is GridMode.RECURSIVE:
        result = recursive_boundary_search(ring, planes, npoints, amplitudes,
                                           nturns=nturns, refpts=refpts, dp=dp,
                                           offset=offset, bounds=bounds,
                                           use_mp=use_mp, verbose=verbose)
    else:
        result = grid_boundary_search(ring, planes, npoints, amplitudes,
                                      nturns=nturns, refpts=refpts, dp=dp,
                                      offset=offset, bounds=bounds,
                                      grid_mode=grid_mode, use_mp=use_mp,
                                      verbose=verbose)
    return result
