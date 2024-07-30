from functools import partial
from .generic import baseplot
from ..lattice import Lattice, Refpts, axis_, plane_
from ..latticetools import LocalOpticsObservable, TrajectoryObservable
from ..latticetools import GeometryObservable, ObservableList


def subscript(plane):
    return "".join(str(axis_(i, "index")) for i in plane)


_pinfo = dict(
    alpha=(r"$\{param}a_{{{plane}}}$", r"$\{param}$", partial(plane_, key="label")),
    beta=(r"$\{param}_{{{plane}}}$", r"$\{param}$ [m]", partial(plane_, key="label")),
    gamma=(r"$\{param}_{{{plane}}}$", r"$\{param}$", partial(plane_, key="label")),
    mu=(r"$\{param}_{{{plane}}}$", r"$\{param}$ [rad]", partial(plane_, key="label")),
    closed_orbit=(r"${{{plane}}}_{{co}}$", "{plane}", partial(axis_, key="code")),
    dispersion=(r"$\eta_{{{plane}}}$", r"$\eta$ [m]", partial(axis_, key="code")),
    s_pos=("s", "s [m]", partial(plane_, key="label")),
    M=(r"${param}_{{{plane}}}$", "{param}", subscript),
    A=(r"${param}_{{{plane}}}$", "{param}", subscript),
    B=(r"${param}_{{{plane}}}$", "{param}", subscript),
    C=(r"${param}_{{{plane}}}$", "{param}", subscript),
    R=(r"${param}_{{{plane}}}$", "{param}", subscript),
    W=(r"${param}_{{{plane}}}$", "{param}", partial(plane_, key="label")),
    Wp=(r"${param}_{{{plane}}}$", "{param}", partial(plane_, key="label")),
    dalpha=(
        r"$\partial \alpha_{{{plane}}}/ \partial \delta$",
        r"$\partial \alpha / \partial \delta$",
        partial(plane_, key="label"),
    ),
    dbeta=(
        r"$\partial \beta_{{{plane}}}/ \partial \delta$",
        r"$\partial \beta/ \partial \delta$ [m]",
        partial(plane_, key="label"),
    ),
    dmu=(
        r"$\partial \mu_{{{plane}}}/ \partial \delta$",
        r"$\partial \mu/ \partial \delta$ [rad]",
        partial(plane_, key="label"),
    ),
    ddispersion=(
        r"$\partial \eta_{{{plane}}}/ \partial \delta$",
        r"$\partial \eta/ \partial \delta$ [m]",
        partial(axis_, key="code"),
    ),
    dR=(
        r"$\partial R_{{{plane}}}/ \partial \delta$",
        r"$\partial R/ \partial \delta$",
        subscript,
    ),
    Trajectory=("{plane}", "{plane}", partial(axis_, key="label")),
    Geometry=("{plane}", "{plane}", partial(plane_, key="label")),
)

_default = ("{param}[{plane}]", "{param}", lambda x: x)


def _plotlabel(param, plane):
    fmt, _, code = _pinfo.get(param, _default)
    return fmt.format(param=param, plane=code(plane))


def _axislabel(param, plane):
    _, fmt, code = _pinfo.get(param, _default)
    return fmt.format(param=param, plane=code(plane))


def plot_optics(ring: Lattice, left: tuple, right: tuple = (), **kwargs):
    """Plot the value of any LocalOpticsObservable

    Parameters:
        ring:   Lattice description.
        left:   Definition of the observables to plot on the left axis
          left[0]: parameter name or user-defined function. See
          :py:class:`.LocalOpticsObservable`
          left[1]: Index in the parameter array or None for scalar values
        right:  Definition of the observables to plot on the right axis

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

    Returns:
        left_axes (Axes):   Main (left) axes
        right_axes (Axes):  Secondary (right) axes or :py:obj:`None`
        synopt_axes (Axes): Synoptic axes

    """

    def plot_function(ring: Lattice, refpts: Refpts, **kwargs):
        def getobs(refpts, param, plane):
            if param == "Trajectory":
                titles.add("Trajectory")
                return TrajectoryObservable(refpts, plane)
            elif param == "Geometry":
                titles.add("Geometry")
                return GeometryObservable(refpts, plane_(plane, "code"))
            else:
                titles.add("Optical functions ")
                return LocalOpticsObservable(refpts, param, plane)

        def process(axes, posdata, data):
            def scan(param, plane):
                axlabs.add(_axislabel(param, plane))
                return _plotlabel(param, plane)

            axlabs = set()
            legend = [scan(*plot) for plot in axes]
            axlabel = ", ".join(axlab for axlab in axlabs)
            return axlabel, posdata, data, legend

        titles = set()
        obs = ObservableList(getobs(refpts, *plot) for plot in left)
        obs.extend(getobs(refpts, *plot) for plot in right)
        obs.append(LocalOpticsObservable(refpts, "s_pos"))
        obs.evaluate(ring, **kwargs)
        leftax = process(left, obs.values[-1], obs.values[:nl])
        if nr > 0:
            rightax = process(right, obs.values[-1], obs.values[nl : nl + nr])
            return ", ".join(t for t in titles), leftax, rightax
        else:
            return ", ".join(t for t in titles), leftax

    nl = len(left)
    nr = len(right)
    return baseplot(ring, plot_function, **kwargs)
