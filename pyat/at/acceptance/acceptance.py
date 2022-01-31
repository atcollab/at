import at
import numpy
from .boundary import GridMode
from .boundary import boundary_search
import multiprocessing
from ..lattice import Lattice


__all__ = ['get_acceptance', 'get_1d_acceptance', 'get_horizontal_acceptance',
           'get_vertical_acceptance', 'get_momentum_acceptance']


def get_acceptance(ring, planes, npoints, amplitudes, nturns=1024,
                   refpts=None, dp=None, offset=None, bounds=None,
                   grid_mode=GridMode.RADIAL, use_mp=False, verbose=True):
    """
    Computes the acceptance at repfts observation points
    Grid Coordiantes ordering is as follows: CARTESIAN: (x,y), RADIAL/RECURSIVE
    (r, theta). Scalar inputs can be used for 1D grid

    PARAMETERS
        ring            ring use for tracking
        planes          max. dimension 2, defines the plane where to search
                        for the acceptance, allowed values are: x,xp,y,yp,dp,ct
        npoints         number of points in each dimension shape (len(planes),)
        amplitudes      max. amplitude  or initial step in RECURSIVE in each
                        dimension
                        shape (len(planes),), for RADIAL/RECURSIVE grid:
                        r = sqrt(x**2+y**2)
        grid_mode       at.GridMode.CARTESIAN: (x,y) grid
                        at.GridMode.RADIAL: (r,theta) grid
                        at.GridMode.RECURSIVE: (r,theta) recursive boundary
                        search

    KEYWORDS
        nturns=1024     Number of turns for the tracking
        refpts=None     Observation refpts, default start of the machine
        dp=None         static momentum offset
        offset=None     initial orbit, default closed orbit
        bounds=None     Allows to define boundaries for the grid default
                        values are:
                        GridMode.CARTESIAN: ((-1,1),(0,1))
                        GridMode.RADIAL/RECURSIVE: ((0,1),(pi,0))
        grid_mode       mode for the gird default GridMode.RADIAL
        use_mp=False    Use python multiprocessing (patpass, default use
                        lattice_pass). In case multi-processing is not
                        enabled GridMode is forced to
                        RECURSIVE (most efficient in single core)
        verbose=True    Print out some inform

    OUTPUT
        Returns 3 lists containing the 2D acceptance, the grid that was
        tracked and the particles of the grid that survived. The length
        of the lists=refpts. In case len(refpts)=1 the acceptance, grid,
        survived arrays are returned directly.
    """
    if not use_mp:
        grid_mode = GridMode.RECURSIVE

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
            cond = grid_mode is GridMode.RADIAL or grid_mode is GridMode.CARTESIAN
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
                                  nturns=nturns, refpts=r, dp=dp,
                                  offset=offset, bounds=bounds,
                                  grid_mode=grid_mode, use_mp=use_mp,
                                  verbose=verbose)
        boundary.append(b)
        survived.append(s)
        grid.append(g)
    if len(rp) == 1:
        return boundary[0], survived[0], grid[0]
    else:
        return boundary, survived, grid


def get_1d_acceptance(ring, plane, resolution, amplitude, nturns=1024, dp=None,
                      refpts=None, grid_mode=GridMode.RADIAL, use_mp=False,
                      verbose=False):
    """
    Computes the 1D acceptance at refpts observation points
    Scalar parameters required

    PARAMETERS
        ring            ring use for tracking
        plane           max. dimension 2, defines the plane where to search
                        for the acceptance, allowed values are: x,xp,y,yp,dp,ct
        resolution      minimum distance between 2 grid points
        amplitude       max. amplitude of the grid or initial step in RECURSIVE

    KEYWORDS
        nturns=1024     Number of turns for the tracking
        dp=0            static momentum offset
        refpts=None     Observation refpts, default start of the machine
        grid_mode       at.GridMode.CARTESIAN/RADIAL: track full vector (default)
                        at.GridMode.RECURSIVE: recursive search
        use_mp=False    Use python multiprocessing (patpass, default use
                        lattice_pass).
                        In case multi-processing is not enabled GridMode is
                        forced to RECURSIVE (most efficient in single core)
        verbose=False   Print out some information


    OUTPUT
        Returns 3 lists containing the 1D acceptance, the grid that was tracked
        and the particles of the grid that survived.
        The length of the lists=refpts. In case len(refpts)=1 the acceptance,
        grid, suvived arrays are returned.
        The boundary output is squeezed to an array with shape (len(refpts),2)
    """
    assert len(numpy.atleast_1d(plane)) == 1, \
        '1D acceptance: single plane required'
    assert numpy.isscalar(resolution), '1D acceptance: scalar args required'
    assert numpy.isscalar(amplitude), '1D acceptance: scalar args required'
    npoint = numpy.ceil(amplitude/resolution)
    b, s, g = get_acceptance(ring, plane, npoint, amplitude,
                             nturns=nturns, dp=dp, refpts=refpts,
                             grid_mode=grid_mode, use_mp=use_mp,
                             verbose=verbose)
    return numpy.squeeze(b), s, g


def get_horizontal_acceptance(ring, *args, **kwargs):
    """
    Computes the 1D horizontal acceptance at refpts observation points
    Scalar parameters required

    PARAMETERS
        ring            ring use for tracking
        resolution      minimum distance between 2 grid points
        amplitude       max. amplitude of the grid or initial step in RECURSIVE

    KEYWORDS
        nturns=1024     Number of turns for the tracking
        dp=0            static momentum offset
        refpts=None     Observation refpts, default start of the machine
        grid_mode       at.GridMode.CARTESIAN/RADIAL: track full vector (default)
                        at.GridMode.RECURSIVE: recursive search
        use_mp=False    Use python multiprocessing (patpass, default use
                        lattice_pass).
                        In case multi-processing is not enabled GridMode is
                        forced to RECURSIVE (most efficient in single core)
        verbose=False   Print out some information


    OUTPUT
        Returns 3 lists containing the 1D acceptance, the grid that was tracked
        and the particles of the grid that survived.
        The length of the lists=refpts. In case len(refpts)=1 the acceptance,
        grid, suvived arrays are returned.
        The boundary output is squeezed to an array with shape (len(refpts),2)
    """
    return get_1d_acceptance(ring, 'x', *args, **kwargs)


def get_vertical_acceptance(ring, *args, **kwargs):
    """
    Computes the 1D vertical acceptance at refpts observation points
    Scalar parameters required

    PARAMETERS
        ring            ring use for tracking
        resolution      minimum distance between 2 grid points
        amplitude       max. amplitude of the grid or initial step in RECURSIVE

    KEYWORDS
        nturns=1024     Number of turns for the tracking
        dp=0            static momentum offset
        refpts=None     Observation refpts, default start of the machine
        grid_mode       at.GridMode.CARTESIAN/RADIAL: track full vector (default)
                        at.GridMode.RECURSIVE: recursive search
        use_mp=False    Use python multiprocessing (patpass, default use
                        lattice_pass).
                        In case multi-processing is not enabled GridMode is
                        forced to RECURSIVE (most efficient in single core)
        verbose=False   Print out some information


    OUTPUT
        Returns 3 lists containing the 1D acceptance, the grid that was tracked
        and the particles of the grid that survived.
        The length of the lists=refpts. In case len(refpts)=1 the acceptance,
        grid, suvived arrays are returned.
        The boundary output is squeezed to an array with shape (len(refpts),2)
    """
    return get_1d_acceptance(ring, 'y', *args, **kwargs)


def get_momentum_acceptance(ring, *args, **kwargs):
    """
    Computes the 1D momentum acceptance at refpts observation points
    Scalar parameters required

    PARAMETERS
        ring            ring use for tracking
        resolution      minimum distance between 2 grid points
        amplitude       max. amplitude of the grid or initial step in RECURSIVE

    KEYWORDS
        nturns=1024     Number of turns for the tracking
        dp=0            static momentum offset
        refpts=None     Observation refpts, default start of the machine
        grid_mode       at.GridMode.CARTESIAN/RADIAL: track full vector (default)
                        at.GridMode.RECURSIVE: recursive search
        use_mp=False    Use python multiprocessing (patpass, default use
                        lattice_pass).
                        In case multi-processing is not enabled GridMode is
                        forced to RECURSIVE (most efficient in single core)
        verbose=False   Print out some information


    OUTPUT
        Returns 3 lists containing the 1D acceptance, the grid that was tracked
        and the particles of the grid that survived.
        The length of the lists=refpts. In case len(refpts)=1 the acceptance,
        grid, suvived arrays are returned.
        The boundary output is squeezed to an array with shape (len(refpts),2)
    """
    return get_1d_acceptance(ring, 'dp', *args, **kwargs)


Lattice.get_acceptance = get_acceptance
Lattice.get_horizontal_acceptance = get_horizontal_acceptance
Lattice.get_vertical_acceptance = get_vertical_acceptance
Lattice.get_momentum_acceptance = get_momentum_acceptance
