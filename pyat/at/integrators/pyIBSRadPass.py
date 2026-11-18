"""
IBS passmethod for Acclerator toolbox(AT) inspired from mbtrack2
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
from at.collective import ibs


def kick(
    bunch, n_bin, sigma_p, sigma_px, sigma_py, revolution_frequency, T_x, T_y, T_p
):
    """
    Applies IBS-induced momentum kicks to the bunch every turn, scaled by growth rates,
    with optional binning to capture non-uniform bunch profiles.
    """

    if n_bin > 1:

        data = bunch[5, :] / clight
        bin_min = np.min(data)
        bin_min = min(bin_min * 0.99, bin_min * 1.01)
        bin_max = np.max(data)
        bin_max = max(bin_max * 0.99, bin_max * 1.01)

        bins = np.linspace(bin_min, bin_max, n_bin + 1)
        center = (bins[1:] + bins[:-1]) / 2
        sorted_index = np.searchsorted(bins, data, side="left")
        sorted_index -= 1
        profile = np.bincount(sorted_index, minlength=n_bin)

        normalized_profile = profile / max(profile)
        Rho = normalized_profile[sorted_index]

    else:

        Rho = 1.0

    N_mp = len(bunch[0, :])
    Delta_pz = (
        sigma_p
        * np.sqrt(np.sqrt(2) * T_p * (1 / revolution_frequency) * Rho)
        * np.random.normal(size=N_mp)
    )
    Delta_px = (
        sigma_px
        * np.sqrt(np.sqrt(2) * T_x * (1 / revolution_frequency) * Rho)
        * np.random.normal(size=N_mp)
    )
    Delta_py = (
        sigma_py
        * np.sqrt(np.sqrt(2) * T_y * (1 / revolution_frequency) * Rho)
        * np.random.normal(size=N_mp)
    )

    bunch[1, :] += Delta_px
    bunch[3, :] += Delta_py
    bunch[4, :] += Delta_pz


def trackFunction(rin, elem=None):

    model = elem.model
    revolution_frequency = elem.revolution_frequency
    n_points = int(elem.n_points)
    n_bin = elem.n_bin
    gamma = elem.gamma
    beta = elem.beta
    r_0 = elem.r_0
    update_turns = int(elem.update_turns)

    ring = elem.ring
    beta_x = elem.beta_x
    beta_y = elem.beta_y
    alphaX = elem.alphaX
    alphaY = elem.alphaY
    dispX = elem.dispX
    disppX = elem.disppX
    dispY = elem.dispY
    disppY = elem.disppY

    dispX_sq_over_beta_x = elem.dispX_sq_over_beta_x
    dispY_sq_over_beta_y = elem.dispY_sq_over_beta_y

    H_x = elem.H_x
    H_y = elem.H_y

    N = ring.beam_current / (qe * revolution_frequency)

    sigma_py, sigma_px, sigma_p, sigma_s, d = ibs.initialize(bunch=rin, ring=ring)

    if elem.Turn_counter % update_turns == 0:

        a, b, q, sigma_H, emitx, emity = ibs.paramerters(
            bunch=rin,
            ring=ring,
            model=model,
            beta_x=beta_x,
            beta_y=beta_y,
            sigma_p=sigma_p,
            dispX_sq_over_beta_x=dispX_sq_over_beta_x,
            dispY_sq_over_beta_y=dispY_sq_over_beta_y,
            H_x=H_x,
            H_y=H_y,
            gamma=gamma,
            beta=beta,
            d=d,
            r_0=r_0,
        )

        if model == "Bane":
            C_log = np.log(q**2 / a**2)
            C_a = a / b

        elif model in ["PM", "PS", "CIMP"]:
            A = (r_0**2 * N) / (
                64 * np.pi**2 * beta**3 * gamma**4 * emitx * emity * sigma_s * sigma_p
            )

        if model in ["PM", "PS"]:

            vabq, v1aq, v1bq = ibs.scatter(
                model=model, n_points=n_points, b=b, a=a, q=q
            )

            T_x, T_y, T_p = ibs.get_scatter_T(
                vabq=vabq,
                v1aq=v1aq,
                v1bq=v1bq,
                A=A,
                model=model,
                sigma_H=sigma_H,
                dispX=dispX,
                dispY=dispY,
                beta_x=beta_x,
                beta_y=beta_y,
                r_0=r_0,
                N=N,
                gamma=gamma,
                beta=beta,
                emitx=emitx,
                emity=emity,
                sigma_s=sigma_s,
                sigma_p=sigma_p,
                H_x=H_x,
                H_y=H_y,
                b=b,
                a=a,
                q=q,
                dispX_sq_over_beta_x=dispX_sq_over_beta_x,
                dispY_sq_over_beta_y=dispY_sq_over_beta_y,
            )

        elif model == "Bane":

            gval = ibs.scatter(model=model, n_points=n_points, b=b, a=a, q=q, C_a=C_a)

            T_x, T_y, T_p = ibs.get_scatter_T(
                gval=gval,
                model=model,
                sigma_H=sigma_H,
                dispX=dispX,
                dispY=dispY,
                beta_x=beta_x,
                beta_y=beta_y,
                r_0=r_0,
                N=N,
                C_log=C_log,
                gamma=gamma,
                beta=beta,
                emitx=emitx,
                emity=emity,
                sigma_s=sigma_s,
                sigma_p=sigma_p,
                H_x=H_x,
                H_y=H_y,
                b=b,
                a=a,
                q=q,
                dispX_sq_over_beta_x=dispX_sq_over_beta_x,
                dispY_sq_over_beta_y=dispY_sq_over_beta_y,
            )

        elif model == "CIMP":

            g_ab, g_ba = ibs.scatter(model=model, n_points=n_points, b=b, a=a, q=q)

            T_x, T_y, T_p = ibs.get_scatter_T(
                g_ab=g_ab,
                g_ba=g_ba,
                A=A,
                model=model,
                sigma_H=sigma_H,
                dispX=dispX,
                dispY=dispY,
                beta_x=beta_x,
                beta_y=beta_y,
                r_0=r_0,
                N=N,
                gamma=gamma,
                beta=beta,
                emitx=emitx,
                emity=emity,
                sigma_s=sigma_s,
                sigma_p=sigma_p,
                H_x=H_x,
                H_y=H_y,
                b=b,
                a=a,
                q=q,
                dispY_sq_over_beta_y=dispY_sq_over_beta_y,
            )

        elem.T_x = T_x
        elem.T_y = T_y
        elem.T_p = T_p

    elem.Turn_counter += 1

    kick(
        rin,
        n_bin,
        sigma_p,
        sigma_px,
        sigma_py,
        revolution_frequency,
        elem.T_x,
        elem.T_y,
        elem.T_p,
    )
