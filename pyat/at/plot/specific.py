"""AT plotting functions"""
from ..lattice import Lattice, Refpts
from ..tracking import internal_lpass
from ..physics import get_optics, Orbit
from .generic import baseplot

# --------- Example 1 --------

# The specific data generating functions do not depend on any graphics-related
# function. So they could be located in any module, for instance
# in lattice/linear.py where the get_optics function is.


# data generating function
def pldata_beta_disp(ring: Lattice, refpts: Refpts, **kwargs):
    """Generates data for plotting beta functions and dispersion"""

    # compute linear optics at the required locations
    data0, _, data = get_optics(ring, refpts=refpts, **kwargs)

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
def plot_beta(ring: Lattice, **kwargs):
    """Plot beta functions and dispersion

    Parameters:
        ring:               Lattice description.

    Keyword Args:
        dp (float):         Momentum deviation.
        dct (float):        Path lengthening. If specified, ``dp`` is
          ignored and the off-momentum is deduced from the path lengthening.
        orbit (Orbit):      Avoids looking for the closed orbit if is
          already known ((6,) array)
        method (Callable):  Method for linear optics (see
          :py:func:`.get_optics`):

          :py:obj:`~.linear.linopt2`: no longitudinal motion, no H/V coupling,

          :py:obj:`~.linear.linopt4`: no longitudinal motion, Sagan/Rubin
          4D-analysis of coupled motion,

          :py:obj:`~.linear.linopt6` (default): with or without longitudinal
          motion, normal mode analysis
        keep_lattice (bool):    Assume no lattice change since the
          previous tracking. Defaults to :py:obj:`False`
        XYStep (float):     Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):     Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`
        twiss_in:           Initial conditions for transfer line optics. Record
          array as output by :py:func:`.linopt2`, :py:func:`.linopt6`, or
          dictionary.
        s_range:            Lattice range of interest, default: unchanged,
          initially set to the full cell.
        axes (tuple[Axes, Optional[Axes]): :py:class:`~matplotlib.axes.Axes`
          for plotting as (primary_axes, secondary_axes).
          Default: create new axes
        slices (int):       Number of slices. Default: 400
        legend (bool):      Show a legend on the plot
        labels (Refpts):    display the name of selected elements.
          Default: :py:obj:`None`
        block (bool):       If :py:obj:`True`, block until the figure is closed.
          Default: :py:obj:`False`
        dipole (dict):      Dictionary of properties overloading the default
          properties of dipole representation. See :py:func:`.plot_synopt`
          for details
        quadrupole (dict):  Same definition as for dipole
        sextupole (dict):   Same definition as for dipole
        multipole (dict):   Same definition as for dipole
        monitor (dict):     Same definition as for dipole
   """
    return baseplot(ring, pldata_beta_disp, **kwargs)

# --------- Example 2 --------


# data generating function
def pldata_linear(ring: Lattice, refpts: Refpts, *keys, **kwargs):
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
def plot_linear(ring: Lattice, *keys, **kwargs):
    """Plot linear optical functions returned by :py:func:`.get_optics`

    Parameters:
        ring:           Lattice description.

    Keyword Args:
        left:           Left axis description as a tuple:
          (key[, indices[, indices]])

          **key**:        'beta', 'closed_orbit',... See :py:func:`.get_optics`

          **indices**:    integer, sequence of integers, or slice

          The number if sequences of indices is data[key].ndim-1
          The number of indices is the number of curves to plot.
          All sequences must have the same length.

        Examples:

          :code:`('beta', [0, 1])`:         beta_x, beta_z

          :code:`('dispersion', 0)`:        eta_x

          :code:`('closed_orbit', [1, 3])`:    x', z'

          :code:`('m44', 2, 2)`:            T33

          :code:`('m44', [0, 0], [0, 1])`:  T11, T12

          :code:`('m44', 2, slice(4))`:     T31, T32, T33, T34 as a single block

          :code:`('m44', [2,2,2,2], [0,1,2,3])`:    T31, T32, T33, T34

        right:              Right axis (optional). See ``left``
        title (str):        Plot title, defaults to "Linear optics"
        dp (float):         Momentum deviation.
        dct (float):        Path lengthening. If specified, ``dp`` is
          ignored and the off-momentum is deduced from the path lengthening.
        orbit (Orbit):      Avoids looking for the closed orbit if is
          already known ((6,) array)
        method (Callable):  Method for linear optics (see
          :py:func:`.get_optics`):

          :py:obj:`~.linear.linopt2`: no longitudinal motion, no H/V coupling,

          :py:obj:`~.linear.linopt4`: no longitudinal motion, Sagan/Rubin
          4D-analysis of coupled motion,

          :py:obj:`~.linear.linopt6` (default): with or without longitudinal
          motion, normal mode analysis
        keep_lattice (bool):    Assume no lattice change since the
          previous tracking. Defaults to :py:obj:`False`
        XYStep (float):     Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):     Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`
        twiss_in:           Initial conditions for transfer line optics. Record
          array as output by :py:func:`~.linear.linopt2`,
          :py:func:`~.linear.linopt6`, or dictionary.
        s_range:            Lattice range of interest, default: unchanged,
          initially set to the full cell.
        axes (tuple[Axes, Optional[Axes]): :py:class:`~matplotlib.axes.Axes`
          for plotting as (primary_axes, secondary_axes).
          Default: create new axes
        slices (int):       Number of slices. Default: 400
        legend (bool):      Show a legend on the plot
        labels (Refpts):    display the name of selected elements.
          Default: :py:obj:`None`
        block (bool):       If :py:obj:`True`, block until the figure is closed.
          Default: :py:obj:`False`
        dipole (dict):      Dictionary of properties overloading the default
          properties of dipole representation. See :py:func:`.plot_synopt`
          for details
        quadrupole (dict):  Same definition as for dipole
        sextupole (dict):   Same definition as for dipole
        multipole (dict):   Same definition as for dipole
        monitor (dict):     Same definition as for dipole
    """
    return baseplot(ring, pldata_linear, *keys, **kwargs)


# --------- Example 3 --------

# Here the data generating function is embedded in the convenience function
def plot_trajectory(ring: Lattice, r_in, nturns: int = 1, **kwargs):
    """Plot a particle's trajectory

    Parameters:
        ring:           Lattice object
        r_in:           (6,n) array: input coordinates of n particles
        nturns:         Number of turns

    Keyword Args:
        keep_lattice (bool):   Assume no lattice change since the
          previous tracking. Defaults to :py:obj:`False`
        s_range:            Lattice range of interest, default: unchanged,
          initially set to the full cell.
        axes (tuple[Axes, Optional[Axes]): :py:class:`~matplotlib.axes.Axes`
          for plotting as (primary_axes, secondary_axes).
          Default: create new axes
        slices (int):       Number of slices. Default: 400
        legend (bool):      Show a legend on the plot
        block (bool):       If :py:obj:`True`, block until the figure is closed.
          Default: :py:obj:`False`
        dipole (dict):      Dictionary of properties overloading the default
          properties of dipole representation. See :py:func:`.plot_synopt`
          for details
        quadrupole (dict):  Same definition as for dipole
        sextupole (dict):   Same definition as for dipole
        multipole (dict):   Same definition as for dipole
        monitor (dict):     Same definition as for dipole
    """
    # noinspection PyShadowingNames
    def pldata_trajectory(ring, refpts, r_in, nturns=1, **kwargs):
        r_out = internal_lpass(ring, r_in, refpts=refpts, nturns=nturns,
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
