"""
IBS module for Acclerator toolbox(AT) inspired from mbtrack2
from Alexin Gamelin , Vadim Gubaidulin (SOLEIL)
Github : https://gitlab.synchrotron-soleil.fr/PA/collective-effects/mbtrack2

"""

import numpy as np
import at
import sys
import scipy.integrate as quad
from scipy.constants import c, elementary_charge, epsilon_0, m_e
from scipy.special import hyp2f1
from scipy.interpolate import interp1d
from at.constants import qe
from at.constants import e_mass
from at.constants import clight


def initialize(bunch, ring):
    """
    calculates bunch parameters at each turn.
    Parameters
    ---------------
    bunch object that is being tracked.
    """

    d = 4 * np.std(bunch[2, :])
    sigma_s = np.std(bunch[5, :] / clight)  # sigma_s in m
    sigma_p = np.std(bunch[4, :])
    sigma_px = np.std(bunch[1, :])
    sigma_py = np.std(bunch[3, :])

    return (sigma_py, sigma_px, sigma_p, sigma_s, d)


def paramerters(
    bunch,
    ring,
    model,
    sigma_p,
    beta_x,
    beta_y,
    dispX_sq_over_beta_x,
    dispY_sq_over_beta_y,
    H_x,
    H_y,
    gamma,
    beta,
    d,
    r_0,
):
    """
    Calculates emittances and derived parameters (a, b, q, sigma_H) based on the selected Intrabeam Scattering (IBS) model.
    These computations are performed every 500 turns by default, or at user-defined intervals specified by the update_turns variable,
    depending on the lattice configuration.

    """

    sigma = at.sigma_matrix(beam=bunch)
    eigs, _ = np.linalg.eig(sigma @ at.jmat(3))
    emits = np.sort(
        np.abs(eigs)[::2]
    )  # order after sorting emits = [emittance_vertical , emittance_horizontal , emittance_longitudinal]
    emitx = emits[1]
    emity = emits[0]

    if model == "PS":
        h = (
            (1 / sigma_p**2)
            + (dispX_sq_over_beta_x / emitx)
            + (dispY_sq_over_beta_y / emity)
        )

        sigma_h = np.sqrt(1 / h)
        sigma_H = sigma_h

    elif model in ["CIMP", "PM", "Bane"]:
        H = (1 / sigma_p**2) + (H_x / emitx) + (H_y / emity)
        sigma_H = np.sqrt(1 / H)

    a = (sigma_H / gamma) * np.sqrt(beta_x / emitx)
    b = (sigma_H / gamma) * np.sqrt(beta_y / emity)
    q = sigma_H * beta * np.sqrt(2 * d / r_0)

    return (a, b, q, sigma_H, emitx, emity)


def scatter(model, n_points, b, a, q, C_a=None):
    """
    Computes IBS scattering integrals using model-specific formulas.
    This calculation is performed every 500 turns by default or at user-defined intervals specified by the update_turns variable,
    depending on the lattice configuration.

    """

    if model in ["PS", "PM"]:
        vabq = np.zeros(n_points, dtype=np.float64)
        v1aq = np.zeros(n_points, dtype=np.float64)
        v1bq = np.zeros(n_points, dtype=np.float64)

        def scattering(u, x, y, z):
            """
            Eq. (17) in:
            L. R. Evans and B. W. Zotter, Intrabeam Scattering in the SPS.
            https://cds.cern.ch/record/126036
            """
            P2 = x**2 + ((1 - x**2) * u**2)
            Q2 = y**2 + ((1 - y**2) * u**2)
            P = np.sqrt(P2)
            Q = np.sqrt(Q2)
            f_abq = (
                8
                * np.pi
                * (1 - 3 * u**2)
                / (P * Q)
                * (2 * np.log(z / 2 * (1 / P + 1 / Q)) - 0.5777777777)
            )

            return f_abq

        for i in range(n_points):
            el_1aq, err = quad.quad(
                scattering, 0, 1, args=(1 / b[i], a[i] / b[i], q[i] / b[i])
            )
            el_1bq, err = quad.quad(
                scattering, 0, 1, args=(1 / a[i], b[i] / a[i], q[i] / a[i])
            )
            el_abq = -(el_1aq * (1 / b[i] ** 2)) - (el_1bq * (1 / a[i] ** 2))
            vabq[i] = el_abq
            v1aq[i] = el_1aq
            v1bq[i] = el_1bq
        return vabq, v1aq, v1bq

    elif model == "Bane":
        gval = np.zeros(n_points, dtype=np.float64)

        def g_func(u, j, C_a):
            """
            Eq. (12) in [2].

            Parameters
            ----------
            u : float
                integration variable.
            j : int
                index.
            C_a : float
                result of a/b

            Returns
            -------
            g_val : array
                Scattering integral value at a given point.

            """
            g_val = ((2 * np.sqrt(C_a[j])) / np.pi) * (
                1 / (np.sqrt(1 + u**2) * np.sqrt(C_a[j] ** 2 + u**2))
            )
            return g_val

        for j in range(n_points):
            reslt, err = quad.quad(g_func, 0, np.inf, args=(j, C_a))
            gval[j] = reslt
        return gval

    elif model == "CIMP":

        def Puv(u, v, x):
            """
            https://dlmf.nist.gov/14.3
            """
            if x < 1:
                val = ((1 + x) / (1 - x)) ** (u / 2) * hyp2f1(
                    v + 1, -v, 1 - u, (0.5 - (0.5 * x))
                )
            else:
                val = ((1 + x) / (x - 1)) ** (u / 2) * hyp2f1(
                    v + 1, -v, 1 - u, (0.5 - (0.5 * x))
                )
            return val

        def g_func(u):
            """
            Eq. (34) in [4].
            """
            x_arg = (1 + u**2) / (2 * u)
            if u >= 1:
                g_val = np.sqrt(np.pi / u) * (
                    (Puv(0, -0.5, x_arg)) + ((3 / 2) * (Puv(-1, -0.5, x_arg)))
                )
            else:
                g_val = np.sqrt(np.pi / u) * (
                    (Puv(0, -0.5, x_arg)) - ((3 / 2) * (Puv(-1, -0.5, x_arg)))
                )
            return g_val

        g_ab = np.zeros(n_points)
        g_ba = np.zeros(n_points)
        for i in range(n_points):
            val_ab = g_func(a[i] / b[i])
            val_ba = g_func(b[i] / a[i])
            g_ab[i] = val_ab
            g_ba[i] = val_ba
        return g_ab, g_ba


def get_scatter_T(
    vabq=None,
    v1aq=None,
    v1bq=None,
    g_ab=None,
    g_ba=None,
    gval=None,
    model=None,
    A=None,
    sigma_H=None,
    dispX=None,
    dispY=None,
    beta_x=None,
    beta_y=None,
    r_0=None,
    N=None,
    C_log=None,
    gamma=None,
    beta=None,
    emitx=None,
    emity=None,
    sigma_s=None,
    sigma_p=None,
    H_x=None,
    H_y=None,
    a=None,
    b=None,
    q=None,
    dispX_sq_over_beta_x=None,
    dispY_sq_over_beta_y=None,
):
    """
    Calculates IBS growth rates (T_x, T_y, T_p) using model-specific scattering integrals and beam parameters.
    This calculation is performed every 500 turns by default or at user-defined intervals specified by the update_turns variable,
    depending on the lattice configuration.

    """

    if model == "PS":
        T_p = A * (vabq * (sigma_H**2 / sigma_p**2))
        T_x = A * (v1bq + (vabq * ((dispX_sq_over_beta_x * sigma_H**2) / (emitx))))
        T_y = A * (v1aq + (vabq * ((dispY_sq_over_beta_y * sigma_H**2) / (emity))))

    elif model == "PM":

        T_p = A * (vabq * (sigma_H**2 / sigma_p**2))
        T_x = A * (v1bq + (vabq * ((H_x * sigma_H**2) / (emitx))))
        T_y = A * (v1aq + (vabq * ((H_y * sigma_H**2) / (emity))))

    elif model == "Bane":
        T_pp = (r_0**2 * N * C_log * sigma_H * gval * (beta_x * beta_y) ** (-1 / 4)) / (
            16 * gamma**3 * emitx ** (3 / 4) * emity ** (3 / 4) * sigma_s * sigma_p**3
        )
        T_p = np.average(T_pp)
        T_x = (sigma_p**2 * H_x * T_pp) / emitx
        T_y = (sigma_p**2 * H_y * T_pp) / emity

    elif model == "CIMP":
        K_a = ((np.log(q**2 / a**2) * g_ba) / a) + ((np.log(q**2 / b**2) * g_ab) / b)

        T_p = 2 * np.pi ** (3 / 2) * A * ((sigma_H**2 / sigma_p**2) * K_a)
        T_x = (
            2
            * np.pi ** (3 / 2)
            * A
            * ((-a * np.log(q**2 / a**2) * g_ba) + (((H_x * sigma_H**2) / emitx) * K_a))
        )
        T_y = (
            2
            * np.pi ** (3 / 2)
            * A
            * ((-b * np.log(q**2 / b**2) * g_ab) + (((H_y * sigma_H**2) / emity) * K_a))
        )

    T_x = np.average(T_x)
    T_y = np.average(T_y)
    T_p = np.average(T_p)

    if T_p <= 0:
        T_p = 0
    if T_x <= 0:
        T_x = 0
    if T_y <= 0:
        T_y = 0

    return T_x, T_y, T_p
