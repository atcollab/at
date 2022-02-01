"""AT plotting functions"""
from at.lattice import Lattice
import matplotlib.pyplot as plt
import numpy


#Function to compute and plot acceptance
def plot_acceptance(ring, *args, **kwargs):
    """
    Computes and plots the acceptance at repfts observation points
    Grid Coordiantes ordering is as follows: GRID: (x,y), RADIAL/RECURSIVE (r, theta).
    Scalar inputs can be used for 1D grid   
    
    PARAMETERS
        ring            ring use for tracking
        planes          max. dimension 2, defines the plane where to search for
                        the acceptance, allowed values are: x,xp,y,yp,dp,ct
        npoints         number of points in each dimension shape (len(planes),)
        amplitudes      max. amplitude  or initial step in RECURSIVE in each dimension
                        shape (len(planes),), for RADIAL/RECURSIVE grid: r = sqrt(x**2+y**2)
        grid_mode       at.GridMode.GRID: (x,y) grid
                        at.GridMode.RADIAL: (r,theta) grid
                        at.GridMode.RECURSIVE: (r,theta) recursive boundary search
                        

    KEYWORDS
        nturns=1024     Number of turns for the tracking
        refpts=None     Observation refpts, default start of the machine
        dp=0            static momentum offset
        offset=None     initial orbit, default no offset
        bounds=None     Allows to define boundaries for the grid default values are:
                        GridMode.GRID: ((-1,1),(0,1))
                        GridMode.RADIAL/RECURSIVE: ((0,1),(pi,0))  
        grid_mode       mode for the gird default GridMode.RADIAL
        use_mp=False    Use python multiprocessing (patpass, default use lattice_pass).
                        In case multi-processing is not enabled GridMode is forced to
                        RECURSIVE (most efficient in single core)
        verbose=True    Print out some inform
        block=False     block execution until the plot is closed, if block=True the figure
                        will automatically close at the end of the execution, to make it
                        persistent plt.show() can be added


    OUTPUT
        Returns 3 lists containing the 2D acceptance, the grid that was tracked and the
        particles of the grid that survived. The length of the lists=refpts. In case
        len(refpts)=1 the acceptance, grid, suvived arrays are returned directly.
    """

    units = {'x':'[m]', 'xp':'[rad]','y':'[m]', 
             'yp':'[rad]', 'dp':'', 'ct':'[m]'}

    refpts = kwargs.pop('refpts',None)
    block = kwargs.pop('block',False)
    if len(numpy.atleast_1d(refpts))>1:
        print('Multiple refpts provided for acceptance plot {0}'.format(refpts))
        print('Using only the first one {0}'.format(refpts[0]))
        kwargs.update('refpts',refpts[0])
    boundary, survived, grid = ring.get_acceptance(*args, **kwargs)
    planes = args[0]
    plt.figure()
    plt.plot(*grid,'.', label='Tracked particles')
    plt.plot(*survived,'.', label='Survived particles')
    if len(planes)==1:
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


Lattice.plot_acceptance = plot_acceptance
