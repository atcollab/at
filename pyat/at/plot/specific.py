"""AT plotting functions"""

from at.plot import baseplot
from at.lattice import get_s_pos
from at.physics import linopt
from at.tracking import lattice_pass

# --------- Example 1 --------

# The specific data generating function does not depend on any
# graphics related function. So it could be located in any module, for instance
# in lattice/linear.py where the linopt function is.


def pldata_beta_disp(ring, refpts, **kwargs):
    """Generates data for plotting beta functions and dispersion"""
    data0, _, _, data = linopt(ring, refpts=refpts, get_chrom=True, **kwargs)
    s_pos = data['s_pos']
    betax = data['beta'][:, 0]
    betaz = data['beta'][:, 1]
    dispersion = data['dispersion'][:, 0]
    left = (s_pos, [betax, betaz], r'$\beta [m]$', [r'$\beta_x$', r'$\beta_z$'])
    right = (s_pos, [dispersion], 'dispersion [m]', ['dispersion'])
    return 'Optical functions', left, right

# A convenience function makes the call easier


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
    def pldata_trajectory(ring, refpts, r_in, **kwargs):
        r_out = lattice_pass(ring, r_in, nturns=nturns, refpts=refpts, **kwargs)
        s_pos = get_s_pos(ring, refpts)
        particles = range(r_out.shape[1])
        xx = [r_out[0, i, :, :] for i in particles]
        zz = [r_out[2, i, :, :] for i in particles]
        xlabels = [r'$x_{0}$'.format(i) for i in particles]
        zlabels = [r'$z_{0}$'.format(i) for i in particles]
        left = (s_pos, xx+zz, 'position [m]', xlabels+zlabels)
        return 'Trajectory', left, None

    return baseplot(ring, pldata_trajectory, r_in, **kwargs)
