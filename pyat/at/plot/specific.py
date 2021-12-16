"""AT plotting functions"""
from at.lattice import Lattice
from at.tracking import lattice_pass
from at.physics import get_optics
from at.plot import baseplot
import matplotlib.pyplot as plt
import numpy

# --------- Example 1 --------

# The specific data generating functions do not depend on any graphics-related
# function. So they could be located in any module, for instance
# in lattice/linear.py where the get_optics function is.


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
              


# data generating function
def pldata_beta_disp(ring, refpts, **kwargs):
    """Generates data for plotting beta functions and dispersion"""

    # compute linear optics at the required locations
    data0, _, data = get_optics(ring, refpts=refpts, get_chrom=True, **kwargs)

    # Extract the plot data
    s_pos = data['s_pos']
    betax = data['beta'][:, 0]
    betaz = data['beta'][:, 1]
    dispersion = data['dispersion'][:, 0]

    # Left axis definition
    left = (r'$\beta$ [m]', s_pos, [betax, betaz],
            [r'$\beta_x$', r'$\beta_z$'])
    # Right axis definition
    right = ('dispersion [m]', s_pos, [dispersion], ['dispersion'])
    return 'Optical functions', left, right


# Convenience function to make the call simpler
def plot_beta(ring, **kwargs):
    """
    Plot beta functions and dispersion

    PARAMETERS
        ring            Lattice object

    KEYWORDS
        dp=0.0          Ignored if radiation is ON. Momentum deviation.
        dct=None        Ignored if radiation is ON. Path lengthening.
                        If specified, dp is ignored and the off-momentum is
                        deduced from the path lengthening.
        method=linopt6  Method used for the analysis of the transfer matrix.
                        See get_optics.
                        linopt6: default
                        linopt2: faster if no longitudinal motion and
                                 no H/V coupling,
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        ddp=1.0E-8      momentum deviation used for computation of
                        chromaticities and dispersion
        twiss_in=None   Initial conditions for transfer line optics. Record
                        array as output by linopt, or dictionary. Keys:
                        'R' or 'alpha' and 'beta'   (mandatory)
                        'closed_orbit',             (default 0)
                        'dispersion'                (default 0)
                        If present, the attribute 'R' will be used, otherwise
                        the attributes 'alpha' and 'beta' will be used. All
                        other attributes are ignored.
   """
    return baseplot(ring, pldata_beta_disp, **kwargs)

# --------- Example 2 --------


# data generating function
def pldata_linear(ring, refpts, *keys, **kwargs):
    """data extraction function for plotting results of get_optics"""

    class Lind(object):
        """helper class for lindata extraction"""
        lab2 = "xz"
        lab6 = ("x", "x'", "z", "z'", "l", "\\delta")
        id6 = "123456"
        params = dict(
            beta=(r'$\beta$ [m]', r'$\beta_{0}$', lab2),
            closed_orbit=('position [m]', r'${0}$', lab6),
            dispersion=('dispersion [m]', r'$\eta_{0}$', lab6),
            alpha=(r'$\alpha$', r'$\alpha_{0}$', lab2),
            mu=(r'Phase advance', r'$\mu_{0}$', lab2),
            gamma=('Gamma', 'Gamma', None),
            M=('Transfert', r'$T_{{{0},{1}}}$', id6)
        )

        @classmethod
        def extract(cls, lindata, key, *idx):
            def it(v):
                try:
                    return iter(v)
                except TypeError:
                    return iter([v])

            axis_title, fmt, convert = cls.params[key]
            indices = list(zip(*(it(i) for i in idx)))
            print(indices)
            datay = (lindata[key][(slice(None),) + ic] for ic in indices)
            labels = (fmt.format(*(convert[i] for i in ic)) for ic in indices)
            return axis_title, lindata['s_pos'], datay, labels

    title = kwargs.pop('title', 'Linear optics')
    data0, _, data = get_optics(ring, refpts=refpts, get_chrom=True, **kwargs)
    return (title,) + tuple(Lind.extract(data, *key) for key in keys)


# Convenience function to make the call simpler
def plot_linear(ring, *keys, **kwargs):
    """
    axleft, axright = plot_linear(ring, left[, right], **keywords
    Plot linear optical functions returned by get_optics

    PARAMETERS
        ring            Lattice object
        left            Left axis description as a tuple:
                        (key[, indices[, indices]])
                          key:        'beta', 'closed_orbit',...
                          indices:    integer, sequence of integers, or slice
                        The number if sequences of indices is data[key].ndim-1
                        The number of indices is the number of curves to plot.
                        All sequences must have the same length.

            Examples:
              ('beta', [0, 1])              beta_x, beta_z
              ('dispersion', 0)             eta_x
              ('closed_orbit'), [1, 3])     x', z'
              ('m44', 2, 2)                 T33
              ('m44', [0, 0], [0, 1])       T11, T12
              ('m44', 2, slice(4))          T31, T32, T33, T34
                                            as a single block
              ('m44', [2,2,2,2], [0,1,2,3]) T31, T32, T33, T34

        right           Right axis (optional)

    KEYWORDS
        title           Plot title, defaults to "Linear optics"
        dp=0.0          Ignored if radiation is ON. Momentum deviation.
        dct=None        Ignored if radiation is ON. Path lengthening.
                        If specified, dp is ignored and the off-momentum is
                        deduced from the path lengthening.
        method=linopt6  Method used for the analysis of the transfer matrix.
                        See get_optics.
                        linopt6: default
                        linopt2: faster if no longitudinal motion and
                                 no H/V coupling,
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        ddp=1.0E-8      momentum deviation used for computation of
                        chromaticities and dispersion
        twiss_in=None   Initial conditions for transfer line optics. Record
                        array as output by linopt, or dictionary. Keys:
                        'R' or 'alpha' and 'beta'   (mandatory)
                        'closed_orbit',             (default 0)
                        'dispersion'                (default 0)
                        If present, the attribute 'R' will be used, otherwise
                        the attributes 'alpha' and 'beta' will be used. All
                        other attributes are ignored.
    """
    return baseplot(ring, pldata_linear, *keys, **kwargs)


# --------- Example 3 --------

# Here the data generating function is embedded in the convenience function
def plot_trajectory(ring, r_in, nturns=1, **kwargs):
    """
    plot a particle's trajectory

    PARAMETERS
        ring            Lattice object
        r_in            6xN array: input coordinates of N particles
        nturns=1        Number of turns

    KEYWORDS
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
    """
    # noinspection PyShadowingNames
    def pldata_trajectory(ring, refpts, r_in, nturns=1, **kwargs):
        r_out = lattice_pass(ring, r_in, refpts=refpts, nturns=nturns,
                             **kwargs)
        s_pos = ring.get_s_pos(refpts)
        particles = range(r_out.shape[1])
        xx = [r_out[0, i, :, :] for i in particles]
        zz = [r_out[2, i, :, :] for i in particles]
        xlabels = [r'$x_{0}$'.format(i) for i in particles]
        zlabels = [r'$z_{0}$'.format(i) for i in particles]
        left = ('position [m]', s_pos, xx+zz, xlabels+zlabels)
        return 'Trajectory', left

    return baseplot(ring, pldata_trajectory, r_in, nturns=nturns, **kwargs)


Lattice.plot_beta = plot_beta
Lattice.plot_trajectory = plot_trajectory
Lattice.plot_acceptance = plot_acceptance
