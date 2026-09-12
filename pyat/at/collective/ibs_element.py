"""
IBS element for Acclerator toolbox(AT) inspired from mbtrack2
from Alexin Gamelin , Vadim Gubaidulin (SOLEIL)
Github : https://gitlab.synchrotron-soleil.fr/PA/collective-effects/mbtrack2

"""

import numpy as np
import at
from enum import IntEnum
from ..lattice import Lattice, AtWarning
from at.lattice import Collective
from at.lattice.elements import _array
from at.lattice.utils import Refpts, uint32_refpts, make_copy
from at.physics import get_optics, avlinopt
from at.constants import clight, qe, e_mass
from typing import Sequence, Optional, Union
import warnings

import scipy.integrate as quad
from scipy.constants import c, elementary_charge, epsilon_0, m_e
from scipy.special import hyp2f1
from scipy.interpolate import interp1d
from at import Element


class IBSElement(Collective, Element):
    """Class to generate an IBS element, inherits from Element"""

    default_pass = {False: "DriftPass", True: "DriftPass"}
    # _conversions = dict(RFCavity._conversions,
    #                    )

    def __init__(self, family_name: str, length: float, ring: Lattice, **kwargs):
        """
        Parameters:
        ring (Lattice):          Lattice object.

        model (str, optional):
           Model to compute the growth rates. Options are "PS", "PM", "Bane", or "CIMP".
           Default is "CIMP".

        n_points (int, optional):
           Number of points at which optics is computed.
           Default is 10000.

        n_bin (int, optional):
           Number of bins for capturing non-uniform beam profiles.
           Default is 100.

        get_opt (str, optional):
           Method to calculate the optics around the lattice object.
           Options: "mbtrack2_eqv", "markers", or "average".
           Default is "markers".

           - "mbtrack2_eqv": Similar to mbtrack2 implementation, uses interpolation to get optics at equidistant points.
           - "markers": Places markers at equidistant points and calculates optics at their positions using get_optics.
           - "average": Places markers at equidistant points and averages optics using avlinopt.

        update_turns (int, optional):
           Number of turns after which growth rates are recomputed.
           Default is 500.
        """

        self.Turn_counter = 0
        self.ring = ring
        self.revolution_frequency = ring.revolution_frequency
        self.r_0 = (qe**2) / (4 * np.pi * epsilon_0 * c**2 * m_e)
        self.gamma = (ring.energy * qe) / (c**2 * m_e)
        self.beta = np.sqrt(1 - (e_mass**2 / ring.energy**2))
        self.T_x = 0
        self.T_y = 0
        self.T_p = 0

        kwargs.setdefault("PassMethod", self.default_pass[True])

        self.model = kwargs.pop("model", "CIMP")
        self.get_opt = kwargs.pop("get_opt", "markers")
        self.n_points = kwargs.pop(
            "n_points", 10000 if self.model in ("mbtrack2_eqv", "markers") else 40000
        )
        self.update_turns = kwargs.pop("update_turns", 500)
        self.n_bin = kwargs.pop("n_b", 100)
        self.circumference = ring.get_s_pos(refpts=len(ring))[0]
        super().__init__(family_name, **kwargs)

    def compute_optics_params(self):

        if self.get_opt == "mbtrack2_eqv":  # most close to the one from mbtrack

            self.ring = self.ring.slice(slices=self.n_points)
            s_points = np.linspace(0, self.circumference, self.n_points)
            _, _, ld_6d = self.ring.disable_6d(copy=True).get_optics(
                refpts=np.arange(len(self.ring) + 1)
            )
            s_pos = self.ring.get_s_pos(refpts=np.arange(len(self.ring) + 1))
            bx_fun = interp1d(s_pos, ld_6d.beta[:, 0])
            by_fun = interp1d(s_pos, ld_6d.beta[:, 1])
            self.beta_x = bx_fun(s_points)
            self.beta_y = by_fun(s_points)
            ax_fun = interp1d(s_pos, ld_6d.alpha[:, 0])
            ay_fun = interp1d(s_pos, ld_6d.alpha[:, 1])
            self.alphaX = ax_fun(s_points)
            self.alphaY = ay_fun(s_points)

            dx_fun = interp1d(s_pos, ld_6d.dispersion[:, 0])
            ddx_fun = interp1d(s_pos, ld_6d.dispersion[:, 1])
            dy_fun = interp1d(s_pos, ld_6d.dispersion[:, 2])
            ddy_fun = interp1d(s_pos, ld_6d.dispersion[:, 3])
            self.dispX = dx_fun(s_points)
            self.disppX = ddx_fun(s_points)
            self.dispY = dy_fun(s_points)
            self.disppY = ddy_fun(s_points)

        elif self.get_opt == "markers":

            marker = at.Marker("ibs_marker")
            marker_list = [marker] * (self.n_points)
            s_points = np.linspace(0, self.circumference, self.n_points, endpoint=False)
            self.ring = self.ring.sbreak(s_points, marker_list)
            refpts = self.ring.get_uint32_index("ibs*")
            _, _, ld_6d = self.ring.disable_6d(copy=True).get_optics(refpts=refpts)
            self.beta_x = ld_6d.beta[:, 0]
            self.beta_y = ld_6d.beta[:, 1]
            self.alphaX = ld_6d.alpha[:, 0]
            self.alphaY = ld_6d.alpha[:, 1]
            self.dispX = ld_6d.dispersion[:, 0]
            self.disppX = ld_6d.dispersion[:, 1]
            self.dispY = ld_6d.dispersion[:, 2]
            self.disppY = ld_6d.dispersion[:, 3]

        elif self.get_opt == "average":

            marker = at.Marker("ibs_marker")
            marker_list = [marker] * (self.n_points)
            s_points = np.linspace(0, self.circumference, self.n_points, endpoint=False)
            self.ring = self.ring.sbreak(s_points, marker_list)
            refpts = self.ring.get_uint32_index("ibs*")
            elemdata, avebeta, _, avedisp, avespos, tune, chrom = self.ring.disable_6d(
                copy=True
            ).avlinopt(dp=0.0, refpts=refpts + 1)
            self.beta_x = avebeta[:, 0]
            self.beta_y = avebeta[:, 1]
            self.alphaX = [e.alpha[0] for e in elemdata]
            self.alphaY = [e.alpha[1] for e in elemdata]
            self.dispX = avedisp[:, 0]
            self.disppX = avedisp[:, 1]
            self.dispY = avedisp[:, 2]
            self.disppY = avedisp[:, 3]

        self.H_x = (1 / self.beta_x) * (
            self.dispX**2
            + ((self.beta_x * self.disppX) + (self.alphaX * self.dispX)) ** 2
        )
        self.H_y = (1 / self.beta_y) * (
            self.dispY**2
            + ((self.beta_y * self.disppY) + (self.alphaY * self.dispY)) ** 2
        )
        self.dispX_sq = self.dispX**2
        self.dispX_sq_over_beta_x = self.dispX_sq / self.beta_x
        self.dispY_sq = self.dispY**2
        self.dispY_sq_over_beta_y = self.dispY_sq / self.beta_y

    def clear_history(self):
        pass
