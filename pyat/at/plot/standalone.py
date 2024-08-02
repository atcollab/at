"""AT plotting functions"""
from __future__ import annotations
from at.lattice import Lattice, axis_
from at.lattice import RFCavity
from at.physics import get_mcf
from at.constants import clight
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy
from numpy import ndarray
from math import sqrt
from fractions import Fraction


# Function to compute and plot acceptance
def plot_acceptance(ring: Lattice, planes, *args, **kwargs):
    # noinspection PyUnresolvedReferences
    r"""Plots the acceptance

    Computes the acceptance at repfts observation points using
    :py:func:`.get_acceptance` and plots the tracked
    and survived particles, and the acceptance boundary.

    Parameters:
        ring:           Lattice definition
        planes:         max. dimension 2, Plane(s) to scan for the acceptance.
          Allowed values are: *'x'*, *'xp'*, *'y'*,
          *'yp'*, *'dp'*, *'ct'*

    Keyword Args:
        acceptance (tuple[ndarray, ndarray, ndarray]): tuple containing
          pre-computed acceptance *(boundary, survived, grid)*
        npoints:        (len(planes),) array: number of points in each
          dimension
        amplitudes:     (len(planes),) array: set the search range:

          * :py:attr:`GridMode.CARTESIAN/RADIAL <.GridMode.RADIAL>`:
            max. amplitude
          * :py:attr:`.GridMode.RECURSIVE`: initial step
        nturns (int):       Number of turns for the tracking
        obspt (Refpts):    Observation points. Default: start of the machine
        dp (float):         Static momentum offset
        offset:             Initial orbit. Default: closed orbit
        bounds:             Defines the tracked range: range=bounds*amplitude.
          It can be used to select quadrants. For example, default values are:

          * :py:attr:`.GridMode.CARTESIAN`: ((-1, 1), (0, 1))
          * :py:attr:`GridMode.RADIAL/RECURSIVE <.GridMode.RADIAL>`: ((0, 1),
            (:math:`\pi`, 0))
        grid_mode (GridMode):   Defines the evaluation grid:

          * :py:attr:`.GridMode.CARTESIAN`: full [:math:`\:x, y\:`] grid
          * :py:attr:`.GridMode.RADIAL`: full [:math:`\:r, \theta\:`] grid
          * :py:attr:`.GridMode.RECURSIVE`: radial recursive search
        use_mp (bool):      Use python multiprocessing (:py:func:`.patpass`,
          default use :py:func:`.lattice_pass`). In case multiprocessing is not
          enabled, *grid_mode* is forced to :py:attr:`.GridMode.RECURSIVE`
          (most efficient in single core)
        verbose (bool):     Print out some information
        divider (int):      Value of the divider used in
          :py:attr:`.GridMode.RECURSIVE` boundary search
        shift_zero:
        start_method (str): Python multiprocessing start method. The default
          :py:obj:`None` uses the python default that is considered safe.
          Available parameters: *'fork'*, *'spawn'*, *'forkserver'*.
          The default for linux is *'fork'*, the default for macOS and
          Windows is *'spawn'*. *'fork'* may be used for macOS to speed up
          the calculation or to solve runtime errors, however  it is
          considered unsafe.

    Returns:
        boundary:   (2,n) array: 2D acceptance
        tracked:    (2,n) array: Coordinates of tracked particles
        survived:   (2,n) array: Coordinates of surviving particles

    Example:
        >>> ring.plot_acceptance(planes, npoints, amplitudes)
        >>> plt.show()
    """
    obspt = kwargs.pop("obspt", None)
    block = kwargs.pop("block", False)
    acceptance = kwargs.pop("acceptance", None)
    if obspt is not None:
        assert numpy.isscalar(obspt), "Scalar value needed for obspt"
    kwargs["refpts"] = obspt
    if acceptance is None:
        boundary, survived, grid = ring.get_acceptance(planes, *args, **kwargs)
    else:
        boundary, survived, grid = acceptance
    plt.figure()
    plt.plot(*grid, ".", label="Tracked particles")
    plt.plot(*survived, ".", label="Survived particles")
    if len(planes) == 1:
        pl0 = axis_(planes[0])
        plt.plot(boundary, numpy.zeros(2), label="Acceptance")
        plt.title("1D {0} acceptance".format(pl0["label"]))
        plt.xlabel("{0}{1}".format(pl0["label"], pl0["unit"]))
    else:
        pl0, pl1 = axis_(planes)
        plt.plot(*boundary, label="Acceptance")
        plt.title("2D {0}-{1} acceptance".format(pl0["label"], pl1["label"]))
        plt.xlabel("{0}{1}".format(pl0["label"], pl0["unit"]))
        plt.xlabel("{0}{1}".format(pl1["label"], pl1["unit"]))
    plt.legend()
    plt.show(block=block)
    return boundary, survived, grid


def plot_geometry(
    ring: Lattice,
    start_coordinates: tuple[float, float, float] = (0, 0, 0),
    centered: bool = False,
    ax: Axes = None,
    **kwargs,
):
    """Compute and plot the 2D ring geometry in cartesian coordinates.

    Parameters:
        ring: Lattice description
        start_coordinates: x,y,angle at starting point
        centered: it True the coordinates origin is the center of the ring
        ax: axes for the plot. If not given, new axes are created

    Keyword arguments are forwarded to the underlying
    :py:func:`~matplotlib.pyplot.plot` function

    Returns:
        geomdata: recarray containing, x, y, angle
        radius: machine radius
        ax: plot axis

    Example:
        >>> ring.plot_geometry()
    """
    if not ax:
        fig, ax = plt.subplots()
    geom, radius = ring.get_geometry(
        start_coordinates=start_coordinates, centered=centered
    )
    ax.plot(
        geom["x"],
        geom["y"],
        "o:",
        linewidth=kwargs.pop("linewidth", 0.5),
        markersize=kwargs.pop("markersize", 2),
        **kwargs,
    )
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_aspect("equal", "box")
    return geom, radius, ax


def plot_sigma(
    sigma,
    axis: tuple[str, str] = ("x", "xp"),
    scale: float = 1.0,
    ax: Axes = None,
    **kwargs,
):
    r"""Plot the projection of the phase space defined by a
    :math:`\Sigma`-matrix on the selected plane.

    Arguments:
        sigma:  :math:`\Sigma`-matrix
        axis:   tuple if indices defining the plane of the :math:`\Sigma`
          projection. Allowed values are: *'x'*, *'xp'*, *'y'*,
          *'yp'*, *'dp'*, *'ct'*. Default: (*'x'*, *'xp'*)
        scale:  Scaling factor for the emittance
        ax: axes for the plot. If not given, new axes are created

    Keyword arguments are forwarded to the underlying
    :py:func:`~matplotlib.pyplot.plot` function
    """
    if not ax:
        fig, ax = plt.subplots()
    ax1, ax2 = axis_(axis)
    axid = axis_(axis, key="index")
    sig22 = sigma[numpy.ix_(axid, axid)]
    eps = sqrt(sig22[0, 0] * sig22[1, 1] - sig22[1, 0] * sig22[0, 1])
    sigx = sqrt(sig22[0, 0])
    tr = numpy.array([[sigx, 0.0], [sig22[0, 1] / sigx, eps / sigx]])
    loop = 2.0 * numpy.pi * numpy.arange(0.0, 1.0, 0.001)
    normcoord = numpy.vstack((numpy.cos(loop), numpy.sin(loop)))
    coord = tr @ normcoord
    line = ax.plot(scale * coord[0, :], scale * coord[1, :], **kwargs)
    ax.set_title("{0}-{1} phase space".format(ax1["label"], ax2["label"]))
    ax.set_xlabel("{0}{1}".format(ax1["label"], ax1["unit"]))
    ax.set_ylabel("{0}{1}".format(ax2["label"], ax2["unit"]))
    return line


def plot_RF_bucket_hamiltonian(
    ring,
    ct_range=None,
    dp_range=None,
    num_points=400,
    num_levels=41,
    plot_separatrix=True,
    **kwargs,
):
    r"""Plot the resulting longitudinal Hamiltonian of a ring (defining the RF
    bucket). The Hamiltonian is calculated by summing all the cavities in the
    ring. Harmonic cavities are supported by the function.

    A perfectly tuned lattice is assumed, the cavities' frequency is nominal
    and the TimeLag is set in a way ensuring ct=0 for the synchronous phase
    by using ring.set_cavity_phase().

    Parameters:
        ring: Lattice description
        ct_range (tuple): Forced lower and upper ct values for the plot.
        Default to :math:`\pm 1.1 \times C / (2h)`
        dp_range (tuple): Forced lower and upper dp values for the plot.
        Default to twice the RF acceptance of the bucket.
        num_points (int): Number of points for 1D grid (ct/dp)
        Default to 400.
        num_levels (int): Number of levels for contour plot. Odd number of
        levels allow to center the colorbar around 0.
        Default to 41.
        plot_separatrix (bool): Flag to plot the separatrix contour
        (:math:`\mathcal{H}=0`).

    Returns:
        CT:   (num_points,num_points) array: ct grid
        DP:    (num_points,num_points) array: dp grid
        hamiltonian:   (num_points,num_points) array: Hamiltonian values
        along (CT,DP) grid
    """
    # Momentum compaction/phase slip factors computed up to third order
    tmp_ring = ring.disable_6d(copy=True)
    alpha = get_mcf(tmp_ring, fit_order=3, n_step=10)

    eta = numpy.zeros(len(alpha))
    eta[0] = alpha[0] - 1 / ring.gamma**2
    eta[1] = 3 * ring.beta**2 / 2 / ring.gamma**2 + alpha[1] - alpha[0] * eta[0]
    eta[2] = (
        -ring.beta**2 * (5 * ring.beta**2 - 1) / (2 * ring.gamma**2)
        + alpha[2]
        - 2 * alpha[0] * alpha[1]
        + alpha[1] / ring.gamma**2
        + alpha[0] ** 2 * eta[0]
        - (3 * ring.beta**2 * alpha[0]) / (2 * ring.gamma**2)
    )

    # (ct, dp) grid calculation (defined around the main RF bucket)
    if ct_range is None:
        ct = numpy.linspace(
            -0.55 * ring.circumference / ring.harmonic_number,
            0.55 * ring.circumference / ring.harmonic_number,
            num=num_points,
        )
    else:
        ct = numpy.linspace(ct_range[0], ct_range[1], num=num_points)
    if dp_range is None:
        U0 = ring.energy_loss
        overvoltage = ring.rf_voltage / U0
        rfa = numpy.sqrt(
            2
            * U0
            / (numpy.pi * alpha[0] * ring.harmonic_number * ring.energy)
            * (numpy.sqrt(overvoltage**2 - 1) - numpy.arccos(1 / overvoltage))
        )
        dp = numpy.linspace(-2 * rfa, 2 * rfa, num=num_points)
    else:
        dp = numpy.linspace(dp_range[0], dp_range[1], num=num_points)
    CT, DP = numpy.meshgrid(ct, dp)

    # Hamiltonian (H=U+T) divided by harmonic number to have
    # U = U(V_rf, h, phi_s)
    # First term of the Hamiltonian
    eta_delta = eta[0] / 2 + eta[1] / 3 * DP + eta[2] / 4 * DP**2
    T = ring.beta**2 * ring.energy * eta_delta * DP**2

    hamiltonian = T
    # Iteration over all lattice cavities
    for rfcav in ring[RFCavity]:
        Voltage = rfcav.Voltage
        HarmNumber = rfcav.HarmNumber
        TimeLag = rfcav.TimeLag

        phi_s = TimeLag * 2 * numpy.pi * rfcav.Frequency / ring.beta / clight
        phi = (
            numpy.pi - phi_s
        ) + CT * 2 * numpy.pi * rfcav.Frequency / ring.beta / clight

        # Second term of the Hamiltonian
        U = (
            Voltage
            / (2 * numpy.pi * HarmNumber)
            * (numpy.cos(phi) - numpy.cos(phi_s) + (phi - phi_s) * numpy.sin(phi_s))
        )
        # Add to total Hamiltonian
        hamiltonian += U

    fig, ax = plt.subplots(1)
    lim_range = numpy.max((numpy.abs(hamiltonian).min(), numpy.abs(hamiltonian).max()))
    levels = numpy.linspace(-lim_range, lim_range, num_levels, endpoint=True)
    co = ax.contourf(CT, DP, hamiltonian, levels, cmap="coolwarm", alpha=0.7)
    # additional contour for visibility
    ax.contour(CT, DP, hamiltonian, levels, cmap="coolwarm")
    if plot_separatrix:
        # separatrix contour
        ax.contour(CT, DP, hamiltonian, [0], colors="black")
        plt.plot([], [], color="black", label="Separatrix")
        ax.legend()
    cb = fig.colorbar(co)
    cb.set_label(r"$\mathcal{H}(ct,\delta)$ [a.u.]", fontsize=18)

    ax.set_xlabel(r"ct [m]")
    ax.set_ylabel(r"$\delta$")

    phi_s = (
        ring.get_rf_timelag()
        * 2
        * numpy.pi
        * ring.get_revolution_frequency()
        * ring.harmonic_number
        / (ring.beta * clight)
    )

    def ct_to_phi(ct):
        return (
            numpy.pi
            - phi_s
            + ct
            / (
                2
                * numpy.pi
                * ring.get_revolution_frequency()
                * ring.harmonic_number
                / clight
            )
        )

    def phi_to_ct(phase):
        return (
            numpy.pi
            - phi_s
            - phase
            * (
                2
                * numpy.pi
                * ring.get_revolution_frequency()
                * ring.harmonic_number
                / clight
            )
        )

    ax2 = ax.secondary_xaxis("top", functions=(phi_to_ct, ct_to_phi))
    ax2.set_xlabel(r"$\phi$ [rad]")

    plt.title(r"$\phi_{RF}$ " + rf"= $\pi -$ {phi_s:.2f}", fontsize=18)

    return CT, DP, hamiltonian


def farey_sequence(nthorder, verbose=False):
    """
    returns the Farey sequency, and the resonance sequence of nth order.

    Arguments:
        nthorder: natural number bigger than 0
    Options:
        verbose: prints extra info. Default: False
    Returns:
        fareyseqfloat: list of elements with the Farey sequence in float format.
            See Eqs.(1,2,3) of [1].
        fareyseqfrac: list of elements with the Farey sequence in frac format.
            See Eqs.(1,2,3) of [1].

    [1] R.Tomas. 'From Farey sequences to resonance diagrams. PRAB 17, 014001 (2014)'
    """
    verboseprint = print if verbose else lambda *a, **k: None
    verboseprint(nthorder)
    farey = []
    fracfarey = []
    af = 0
    bf = 1
    cf = 1
    df = nthorder
    farey.append(0)
    farey.append(1 / df)
    fracfarey.append(Fraction(0))
    fracfarey.append(Fraction(1, df))
    idx = 0
    while (farey[-1] < 1) and (idx < 100):
        idx += 1
        caux = numpy.floor((nthorder + bf) / df) * cf - af
        daux = numpy.floor((nthorder + bf) / df) * df - bf
        af = cf
        bf = df
        cf = int(caux)
        df = int(daux)
        farey.append(cf / df)
        fracfarey.append(Fraction(cf, df))
    verboseprint(f"farey_float{nthorder}= {farey}")
    verboseprint(f"farey_frac{nthorder} = {fracfarey}")
    return farey, fracfarey


def plot_tune2D_resonances(
    orders=[1, 2, 3],
    period=1,
    window=numpy.array([0, 1, 0, 1]),
    verbose=False,
    **kwargs,
):
    """
    This function plot the tune 2D resonances for a given order, period and window.

    # the resonance equation
    # int_the_res[0]*nux + int_the_res[1]*nuy = int_res

    """
    # 2024jul31 oblanco at ALBA CELLS

    # verboseprint to check flag only once
    verboseprint = print if verbose else lambda *a, **k: None
    block = kwargs.pop("block", False)

    listresonancestoplot = list(numpy.array(orders) - 1)
    theperiod = period
    verboseprint(f"The period is {theperiod}")
    verboseprint(f"The window is {window}")

    if sum(n < 0 for n in listresonancestoplot):
        print("Error, list includes negative resonances")

    # check the window
    if window[0] == window[1]:
        print("horizontal coordinates must be different")
        exit()
    if window[2] == window[3]:
        print("vertical coordinates must be different")
        exit()
    if window[1] < window[0]:
        print("Swapping horizontal coordinates")
        window[0], window[1] = window[1], window[0]
    if window[3] < window[2]:
        print("Swapping vertical coordinates")
        window[2], window[3] = window[3], window[2]
    # get xlimits and ylimits
    the_axeslims = window.reshape((2, 2))

    # horizontal and vertical borders
    borders = numpy.eye(2)

    maxreson2calc = numpy.max(listresonancestoplot) + 1
    theorder = maxreson2calc  # ??? unnecessary
    verboseprint(f"Farey max order={theorder}")

    # get the Farey collection
    fareycollectionfloat = []
    fareycollectionfrac = []
    for nthorder in range(1, theorder + 1):
        farey, fracfarey = farey_sequence(nthorder)
        fareycollectionfloat.append(farey.copy())
        fareycollectionfrac.append(fracfarey.copy())
    verboseprint(f"the Farey collection is {fareycollectionfloat}")

    # plot configuration
    fig = plt.figure()
    # window min/max,horizontal and vertical
    plt.xlim(the_axeslims[0, :])
    plt.ylim(the_axeslims[1, :])
    # min/max to plot lines with slopes
    minx = numpy.floor(the_axeslims[0, 0])
    minx = minx - theperiod - numpy.mod(minx, theperiod)
    maxx = numpy.ceil(the_axeslims[0, 1])
    maxx = maxx + theperiod - numpy.mod(maxx, theperiod)
    minmaxxdist = maxx - minx
    miny = numpy.floor(the_axeslims[1, 0])
    miny = miny - theperiod - numpy.mod(miny, theperiod)
    maxy = numpy.ceil(the_axeslims[1, 1])
    maxy = maxy + theperiod - numpy.mod(maxy, theperiod)
    minmaxydist = maxy - miny
    minmaxdist = max([minmaxxdist, minmaxydist])

    # dictionary with line properties
    widthmod = 5  # ??? maybe a variable
    mypalettecolor = {
        0: "k",
        1: "b",
        2: "r",
        3: "g",
        4: "m",
        5: "c",
        6: "y",
        7: "darkcyan",
        8: "lightgreen",
        9: (0.1, 0.1, 0.1),
        10: (0.1, 0.1, 0.1),
        11: (0.2, 0.2, 0.2),
        12: (0.3, 0.3, 0.3),
        13: (0.4, 0.4, 0.4),
        14: (0.5, 0.5, 0.5),
        15: (0.6, 0.6, 0.6),
    }
    mypalettestyle = {0: "-", 1: "--"}
    prop1n = dict(
        color=mypalettecolor[0],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(4, widthmod),
        label="axvline - full height",
    )
    prop2n = dict(
        color=mypalettecolor[1],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(3, widthmod),
        label="axvline - full height",
    )
    prop3n = dict(
        color=mypalettecolor[2],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(2, widthmod),
        label="axvline - full height",
    )
    prop4n = dict(
        color=mypalettecolor[3],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(2, widthmod),
        label="axvline - full height",
    )
    prop5n = dict(
        color=mypalettecolor[4],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop6n = dict(
        color=mypalettecolor[5],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop7n = dict(
        color=mypalettecolor[6],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop8n = dict(
        color=mypalettecolor[7],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop9n = dict(
        color=mypalettecolor[8],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop10n = dict(
        color=mypalettecolor[9],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop11n = dict(
        color=mypalettecolor[10],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop12n = dict(
        color=mypalettecolor[11],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop13n = dict(
        color=mypalettecolor[12],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop14n = dict(
        color=mypalettecolor[13],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop15n = dict(
        color=mypalettecolor[14],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop1s = dict(
        color=mypalettecolor[0],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(4, widthmod),
        label="axvline - full height",
    )
    prop2s = dict(
        color=mypalettecolor[1],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(3, widthmod),
        label="axvline - full height",
    )
    prop3s = dict(
        color=mypalettecolor[2],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(2, widthmod),
        label="axvline - full height",
    )
    prop4s = dict(
        color=mypalettecolor[3],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(2, widthmod),
        label="axvline - full height",
    )
    prop5s = dict(
        color=mypalettecolor[4],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop6s = dict(
        color=mypalettecolor[5],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop7s = dict(
        color=mypalettecolor[6],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop8s = dict(
        color=mypalettecolor[7],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop9s = dict(
        color=mypalettecolor[8],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop10s = dict(
        color=mypalettecolor[9],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop11s = dict(
        color=mypalettecolor[10],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop12s = dict(
        color=mypalettecolor[11],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop13s = dict(
        color=mypalettecolor[12],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop14s = dict(
        color=mypalettecolor[13],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    prop15s = dict(
        color=mypalettecolor[14],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="axvline - full height",
    )
    # assemble default dictionary by normal or skew
    propn = {
        0: prop1n.copy(),
        1: prop2n.copy(),
        2: prop3n.copy(),
        3: prop4n.copy(),
        4: prop5n.copy(),
        5: prop6n.copy(),
        6: prop7n.copy(),
        7: prop8n.copy(),
        8: prop9n.copy(),
        9: prop10n.copy(),
        10: prop11n.copy(),
        11: prop12n.copy(),
        12: prop13n.copy(),
        13: prop14n.copy(),
        14: prop15n.copy(),
    }
    props = {
        0: prop1s.copy(),
        1: prop2s.copy(),
        2: prop3s.copy(),
        3: prop4s.copy(),
        4: prop5s.copy(),
        5: prop6s.copy(),
        6: prop7s.copy(),
        7: prop8s.copy(),
        8: prop9s.copy(),
        9: prop10s.copy(),
        10: prop11s.copy(),
        11: prop12s.copy(),
        12: prop13s.copy(),
        13: prop14s.copy(),
        14: prop15s.copy(),
    }
    # assemble final dict
    lprop = {0: propn.copy(), 1: props.copy()}

    # we only need to points to define a line
    nauxpoints = 2

    # start to check the Farey collection, starting with 0
    collectaux1 = [0]
    for nthorder in range(0, maxreson2calc):
        verboseprint(f"nthorder={nthorder}")
        collectaux2 = fareycollectionfloat[nthorder]
        thesteps = list(set(collectaux2) - set(collectaux1))
        verboseprint(thesteps)
        chosenstep = min(thesteps)
        verboseprint(f"chosenstep={chosenstep}")
        collectaux1 = collectaux2
        if not (nthorder in listresonancestoplot):
            continue
        else:
            # increase step by the period in straight resonances
            verboseprint("enter plotting straight lines")
            for iaux in numpy.arange(minx, maxx + 0.000001, theperiod * chosenstep):
                plt.axvline(x=iaux, **lprop[0][nthorder].copy())
            for iaux in numpy.arange(miny, maxy + 0.000001, theperiod * chosenstep):
                plt.axhline(
                    y=iaux, **lprop[numpy.mod(nthorder + 1, 2)][nthorder].copy()
                )
            for iaux in range(nthorder):
                verboseprint(f"enter plotting diagonals {iaux}")
                aeq = iaux + 1
                beq = nthorder + 1 - aeq
                chosenslope = -aeq / beq
                diagstep = 1 / beq
                verboseprint(f"chosen slope={chosenslope}")
                verboseprint(f"chosen diagstep={diagstep}")
                # get the vertical limits with diagonals
                a1 = (
                    -numpy.ceil(2 * (minmaxxdist + minmaxydist) / (-chosenslope)) + miny
                )
                a2 = numpy.ceil(2 * (minmaxxdist + minmaxydist) / (-chosenslope)) + maxy
                verboseprint(f"minx={minx},maxx={maxx},minmaxxdist={minmaxxdist}")
                verboseprint(f"miny={miny},maxy={maxy},minmaxydist={minmaxydist}")
                verboseprint(f"a1={a1},a2={a2}")
                xaux = numpy.linspace(minx, maxx, nauxpoints)
                verboseprint(f"xaux={xaux}")
                for istep in numpy.arange(0, a2 - a1 + 0.0001, diagstep):
                    y1line = chosenslope * (xaux - minx) + theperiod * istep + miny
                    y2line = -chosenslope * (xaux - minx) - theperiod * istep + maxy
                    verboseprint(f"y1line={y1line},y2line={y2line}")
                    plt.plot(xaux, y1line, **lprop[numpy.mod(beq, 2)][nthorder])
                    plt.plot(xaux, y2line, **lprop[numpy.mod(beq, 2)][nthorder])

    plt.xlabel(r"$\nu_x$")
    plt.ylabel(r"$\nu_y$")
    plt.show(block=block)


Lattice.plot_acceptance = plot_acceptance
Lattice.plot_geometry = plot_geometry
Lattice.plot_RF_bucket_hamiltonian = plot_RF_bucket_hamiltonian
