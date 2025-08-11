"""AT plotting functions"""

from __future__ import annotations

from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from at.constants import clight
from at.lattice import Lattice, axis_
from at.lattice import RFCavity
from at.physics import get_mcf


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
        assert np.isscalar(obspt), "Scalar value needed for obspt"
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
        plt.plot(boundary, np.zeros(2), label="Acceptance")
        plt.title(f"1D {pl0['label']} acceptance")
        plt.xlabel(f"{pl0['label']}{pl0['unit']}")
    else:
        pl0, pl1 = axis_(planes)
        plt.plot(*boundary, label="Acceptance")
        plt.title(f"2D {pl0['label']}-{pl1['label']} acceptance")
        plt.xlabel(f"{pl0['label']}{pl0['unit']}")
        plt.xlabel(f"{pl1['label']}{pl1['unit']}")
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
    sig22 = sigma[np.ix_(axid, axid)]
    eps = sqrt(sig22[0, 0] * sig22[1, 1] - sig22[1, 0] * sig22[0, 1])
    sigx = sqrt(sig22[0, 0])
    tr = np.array([[sigx, 0.0], [sig22[0, 1] / sigx, eps / sigx]])
    loop = 2.0 * np.pi * np.arange(0.0, 1.0, 0.001)
    normcoord = np.vstack((np.cos(loop), np.sin(loop)))
    coord = tr @ normcoord
    line = ax.plot(scale * coord[0, :], scale * coord[1, :], **kwargs)
    ax.set_title(f"{ax1['label']}-{ax2['label']} phase space")
    ax.set_xlabel(f"{ax1['label']}{ax1['unit']}")
    ax.set_ylabel(f"{ax2['label']}{ax2['unit']}")
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
    by using :py:func:`.set_cavity_phase`.

    Parameters:
        ring: Lattice description
        ct_range (tuple): Forced lower and upper ct values for the plot.
            Default to :math:`\pm 1.1 \times C / (2h)`
        dp_range (tuple): Forced lower and upper dp values for the plot.
            Default to twice the RF acceptance of the bucket.
        num_points (int): Number of points for 1D grid (ct/dp)
            Default to 400.
        num_levels (int): Number of levels for contour plot. Odd number of
            levels allow to center the colorbar around 0. Default to 41.
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

    eta = np.zeros(len(alpha))
    eta[0] = alpha[0] - 1 / ring.gamma**2
    eta[1] = 3 * ring.beta**2 / 2 / ring.gamma**2 + alpha[1] - alpha[0] * eta[0]
    eta[2] = (
        -(ring.beta**2) * (5 * ring.beta**2 - 1) / (2 * ring.gamma**2)
        + alpha[2]
        - 2 * alpha[0] * alpha[1]
        + alpha[1] / ring.gamma**2
        + alpha[0] ** 2 * eta[0]
        - (3 * ring.beta**2 * alpha[0]) / (2 * ring.gamma**2)
    )

    # (ct, dp) grid calculation (defined around the main RF bucket)
    if ct_range is None:
        ct = np.linspace(
            -0.55 * ring.circumference / ring.harmonic_number,
            0.55 * ring.circumference / ring.harmonic_number,
            num=num_points,
        )
    else:
        ct = np.linspace(ct_range[0], ct_range[1], num=num_points)
    if dp_range is None:
        U0 = ring.energy_loss
        overvoltage = ring.rf_voltage / U0
        rfa = np.sqrt(
            2
            * U0
            / (np.pi * alpha[0] * ring.harmonic_number * ring.energy)
            * (np.sqrt(overvoltage**2 - 1) - np.arccos(1 / overvoltage))
        )
        dp = np.linspace(-2 * rfa, 2 * rfa, num=num_points)
    else:
        dp = np.linspace(dp_range[0], dp_range[1], num=num_points)
    CT, DP = np.meshgrid(ct, dp)

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

        phi_s = TimeLag * 2 * np.pi * rfcav.Frequency / ring.beta / clight
        phi = (np.pi - phi_s) + CT * 2 * np.pi * rfcav.Frequency / ring.beta / clight

        # Second term of the Hamiltonian
        U = (
            Voltage
            / (2 * np.pi * HarmNumber)
            * (np.cos(phi) - np.cos(phi_s) + (phi - phi_s) * np.sin(phi_s))
        )
        # Add to total Hamiltonian
        hamiltonian += U

    fig, ax = plt.subplots(1)
    lim_range = np.max((np.abs(hamiltonian).min(), np.abs(hamiltonian).max()))
    levels = np.linspace(-lim_range, lim_range, num_levels, endpoint=True)
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
        * np.pi
        * ring.get_revolution_frequency()
        * ring.harmonic_number
        / (ring.beta * clight)
    )

    def ct_to_phi(ct):
        return (
            np.pi
            - phi_s
            + ct
            / (
                2
                * np.pi
                * ring.get_revolution_frequency()
                * ring.harmonic_number
                / clight
            )
        )

    def phi_to_ct(phase):
        return (
            np.pi
            - phi_s
            - phase
            * (
                2
                * np.pi
                * ring.get_revolution_frequency()
                * ring.harmonic_number
                / clight
            )
        )

    ax2 = ax.secondary_xaxis("top", functions=(phi_to_ct, ct_to_phi))
    ax2.set_xlabel(r"$\phi$ [rad]")

    plt.title(r"$\phi_{RF}$ " + rf"= $\pi -$ {phi_s:.2f}", fontsize=18)

    return CT, DP, hamiltonian


Lattice.plot_acceptance = plot_acceptance
Lattice.plot_geometry = plot_geometry
Lattice.plot_RF_bucket_hamiltonian = plot_RF_bucket_hamiltonian
