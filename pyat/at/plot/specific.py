"""AT plotting functions"""
from at.lattice import Lattice
from at.tracking import lattice_pass
from at.physics import linopt
from at.plot import baseplot

# --------- Example 1 --------

# The specific data generating functions do not depend on any graphics-related
# function. So they could be located in any module, for instance
# in lattice/linear.py where the linopt function is.


# data generating function
def pldata_beta_disp(ring, refpts, **kwargs):
    """Generates data for plotting beta functions and dispersion"""

    # compute linear optics at the required locations
    data0, _, _, data = linopt(ring, refpts=refpts, get_chrom=True, **kwargs)

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
        dp=0.0          momentum deviation.
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        ddp=1.0E-8      momentum deviation used for computation of
                        chromaticities and dispersion
        coupled=True    if False, simplify the calculations by assuming
                        no H/V coupling
    """
    return baseplot(ring, pldata_beta_disp, **kwargs)

# --------- Example 2 --------


# data generating function
def pldata_linear(ring, refpts, *keys, **kwargs):
    """data extraction function for plotting results of linopt"""

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
            A=('A', r'$A_{{{0}{1}}}$', id6),
            B=('B', r'$B_{{{0}{1}}}$', id6),
            C=('C', r'$C_{{{0}{1}}}$', id6),
            m44=('Transfert', r'$T_{{{0},{1}}}$', id6)
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
    data0, _, _, data = linopt(ring, refpts=refpts, get_chrom=True, **kwargs)
    return (title,) + tuple(Lind.extract(data, *key) for key in keys)


# Convenience function to make the call simpler
def plot_linear(ring, *keys, **kwargs):
    """
    axleft, axright = plot_linear(ring, left[, right], **keywords
    Plot linear optical functions returned by linopt

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
        dp=0.0          momentum deviation.
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        ddp=1.0E-8      momentum deviation used for computation of
                        chromaticities and dispersion
        coupled=True    if False, simplify the calculations by assuming
                        no H/V coupling
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
