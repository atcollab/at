"""Intra-beam scattering (IBS).

Analytical IBS growth rates, equilibrium emittances with radiation damping,
and a collective element applying the corresponding random momentum kicks
during tracking.

This follows the xsuite implementation (xfields.ibs), for bunched beams:

* ``"Nagaitsev"``: S. Nagaitsev, PRSTAB 8, 064403 (2005). Fast, but
  ignores the vertical dispersion,
* ``"Bjorken-Mtingwa"`` (or ``"B&M"``): J.D. Bjorken, S.K. Mtingwa,
  Part. Accel. 13, 115 (1983), as modified in MAD-X to include the vertical
  dispersion: F. Antoniou, F. Zimmermann, CERN-ATS-2012-066,
* Coulomb logarithm computed as in MAD-X (``twclog``).

The kick model follows R. Bruce et al., PRSTAB 13, 091001 (2010), as also
implemented in mbtrack2, xfields and elegant.
"""

from __future__ import annotations

__all__ = ["ibs_rates", "ibs_equilibrium", "IBSElement"]

import numpy as np
from scipy.constants import epsilon_0, hbar
from scipy.integrate import quad_vec
from scipy.special import elliprd

from ..constants import clight, qe
from ..lattice import All, Collective, Element, Lattice, random

_MODELS = {"nagaitsev": "Nagaitsev", "bjorken-mtingwa": "Bjorken-Mtingwa",
           "b&m": "Bjorken-Mtingwa"}


def _check_model(model: str) -> str:
    try:
        return _MODELS[model.lower()]
    except KeyError:
        raise ValueError(
            f"Unknown IBS model {model!r}, choose 'Nagaitsev' or 'Bjorken-Mtingwa'"
        ) from None


def _ibs_optics(ring: Lattice, step: float) -> np.ndarray:
    """Optics sampled around the ring.

    Returns:
        optics: (9, N) array of integration weights normalised to 1,
          beta_x, beta_y, alpha_x, alpha_y, D_x, D'_x, D_y, D'_y
    """
    lat = ring.disable_6d(copy=True).slice(size=step)
    _, _, ld = lat.get_optics(refpts=All)
    s = ld.s_pos
    w = np.zeros_like(s)  # trapezoidal weights
    w[:-1] += np.diff(s)
    w[1:] += np.diff(s)
    return np.vstack((w / w.sum(), ld.beta.T, ld.alpha.T, ld.dispersion.T))


def _beam_constants(ring: Lattice) -> np.ndarray:
    """Energy [eV], rest energy [eV], charge and classical radius [m]."""
    particle = ring.particle
    if particle.rest_energy == 0.0:
        raise ValueError(
            "IBS needs a massive particle: set ring.particle (e.g. 'electron')"
        )
    r0 = qe * particle.charge**2 / (4.0 * np.pi * epsilon_0 * particle.rest_energy)
    return np.array([ring.energy, particle.rest_energy, particle.charge, r0])


def _coulomb_log(optics, consts, npart, emitx, emity, sigma_e, bunch_length):
    """Coulomb logarithm, as in MAD-X ``twclog``."""
    weights, betx, bety, _, _, dx, _, dy, _ = optics
    energy, mass, charge, _ = consts
    gamma = energy / mass
    bxbar, bybar, dxbar, dybar = weights @ np.stack((betx, bety, dx, dy), axis=1)
    etrans = 5.0e8 * (gamma * energy - mass) * 1.0e-9 * emitx / bxbar
    tempev = 2.0 * etrans
    # Beam sizes in cm
    sigx = 100.0 * np.sqrt(emitx * bxbar + (dxbar * sigma_e) ** 2)
    sigy = 100.0 * np.sqrt(emity * bybar + (dybar * sigma_e) ** 2)
    sigt = 100.0 * bunch_length
    density = npart / (8.0 * np.pi**1.5 * sigx * sigy * sigt)
    debye_length = 743.4 * np.sqrt(tempev / density) / abs(charge)
    rmincl = 1.44e-7 * charge**2 / tempev
    hbar_c = hbar * clight / qe * 1.0e-9  # [GeV.m]
    rminqm = hbar_c * 1.0e5 / (2.0 * np.sqrt(2.0e-3 * etrans * mass * 1.0e-9))
    return np.log(min(sigx, debye_length) / max(rmincl, rminqm))


def _nagaitsev(optics, gamma, emitx, emity, sigma_e):
    """Ring averages of the Nagaitsev integrands, Eqs. (30-35)."""
    weights, betx, bety, alfx, _, dx, dpx, dy, _ = optics
    sigx = np.sqrt(betx * emitx + (dx * sigma_e) ** 2)
    sigy = np.sqrt(bety * emity + (dy * sigma_e) ** 2)
    phix = dpx + alfx * dx / betx
    ax = betx / emitx
    ay = bety / emity
    a_s = ax * (dx**2 / betx**2 + phix**2) + 1.0 / sigma_e**2
    a1 = 0.5 * (ax + gamma**2 * a_s)
    a2 = 0.5 * (ax - gamma**2 * a_s)
    sqrt_term = np.sqrt(a2**2 + gamma**2 * ax**2 * phix**2)
    lambda1 = ay
    lambda2 = a1 + sqrt_term
    lambda3 = a1 - sqrt_term
    r1 = elliprd(1.0 / lambda2, 1.0 / lambda3, 1.0 / lambda1) / lambda1
    r2 = elliprd(1.0 / lambda3, 1.0 / lambda1, 1.0 / lambda2) / lambda2
    r3 = (3.0 * np.sqrt(lambda1 * lambda2 / lambda3)
          - lambda1 * r1 / lambda3 - lambda2 * r2 / lambda3)
    sp = 0.5 * gamma**2 * (2.0 * r1 - r2 * (1.0 - 3.0 * a2 / sqrt_term)
                           - r3 * (1.0 + 3.0 * a2 / sqrt_term))
    sx = 0.5 * (2.0 * r1 - r2 * (1.0 + 3.0 * a2 / sqrt_term)
                - r3 * (1.0 - 3.0 * a2 / sqrt_term))
    sxp = 3.0 * gamma**2 * phix**2 * ax * (r3 - r2) / sqrt_term
    ix = betx * (sx + sp * (dx**2 / betx**2 + phix**2) + sxp)
    iy = bety * (r2 + r3 - 2.0 * r1)
    iz = sp
    return (weights / (sigx * sigy)) @ np.stack((ix, iy, iz), axis=1)


def _bjorken_mtingwa(optics, gamma, emitx, emity, sigma_e):
    """Integrals of Eq. (8) of the MAD-X note, with the factors in brackets.

    The terms of Table 1 are multiplied by the bracket factors, which avoids
    divisions by H_x."""
    weights, betx, bety, alfx, alfy, dx, dpx, dy, dpy = optics
    g2 = gamma**2
    s2 = 1.0 / sigma_e**2
    gx = betx / emitx
    gy = bety / emity
    phix = dpx + alfx * dx / betx
    phiy = dpy + alfy * dy / bety
    hx = (dx**2 + betx**2 * phix**2) / betx / emitx  # H_x / emitx
    hy = (dy**2 + bety**2 * phiy**2) / bety / emity  # H_y / emity
    ry = hy * emity / bety  # H_y / beta_y
    fx = gx**2 * phix**2
    fy = gy**2 * phiy**2
    hsum = hx + hy + s2
    dsum = g2 * (dx**2 / betx / emitx + dy**2 / bety / emity + s2)

    a = g2 * hsum + gx + gy
    b = (gx + gy) * dsum + gx * gy * (g2 * (phix**2 + phiy**2) + 1.0)
    c = gx * gy * dsum
    ax = (g2 * hx * (2.0 * g2 * hsum - 2.0 * gx - gy) - g2 * gx * hy
          + gx * (2.0 * gx - gy - g2 * s2) + 6.0 * g2 * fx)
    bx = (g2 * hx * ((gx + gy) * g2 * hsum - g2 * (fx + fy) + gx * (gx - 4.0 * gy))
          + gx * (g2 * s2 * (gx - 2.0 * gy) + gx * gy * (1.0 + 6.0 * g2 * phix**2)
                  + g2 * (2.0 * fy - fx))
          + g2 * gx * hy * (gx - 2.0 * gy))
    ay = gy * (-g2 * (hx + 2.0 * hy + gx * ry + s2)
               + 2.0 * g2**2 * ry * hsum - (gx - 2.0 * gy) + 6.0 * g2 * gy * phiy**2)
    by = gy * (g2 * (gy - 2.0 * gx) * (hx + s2) + g2 * hy * (gy - 4.0 * gx)
               + gx * gy + g2 * (2.0 * fx - fy)
               + g2**2 * ry * (gx + gy) * hsum - g2**2 * ry * (fx + fy)
               + 6.0 * g2 * phiy**2 * gx * gy)
    az = g2 * s2 * (2.0 * g2 * hsum - gx - gy)
    bz = g2 * s2 * ((gx + gy) * g2 * hsum - 2.0 * gx * gy - g2 * (fx + fy))
    num_a = np.concatenate((ax, ay, az))
    num_b = np.concatenate((bx, by, bz))
    a, b, c = np.tile(a, 3), np.tile(b, 3), np.tile(c, 3)

    def integrand(lam):
        return np.sqrt(lam) * (num_a * lam + num_b) / (
            lam**3 + a * lam**2 + b * lam + c
        ) ** 1.5

    # Integration over [1, 1e16] split in decades, as in MAD-X and xsuite
    result = sum(quad_vec(integrand, 10.0**i, 10.0**(i + 1))[0] for i in range(16))
    return result.reshape(3, -1) @ weights


def _growth_rates(optics: np.ndarray, model: str, consts: np.ndarray,
                  npart: float, emitx: float, emity: float,
                  sigma_e: float, bunch_length: float) -> np.ndarray:
    """Emittance growth rates [1/s] for given beam parameters."""
    energy, mass, _, r0 = consts
    gamma = energy / mass
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    clog = _coulomb_log(optics, consts, npart, emitx, emity, sigma_e, bunch_length)
    if model == "Nagaitsev":
        cst = npart * r0**2 * clight * clog / (
            12.0 * np.pi * beta**3 * gamma**5 * bunch_length
        )
        integrals = _nagaitsev(optics, gamma, emitx, emity, sigma_e)
    else:
        cst = npart * r0**2 * clight * clog / (
            8.0 * np.pi * beta**3 * gamma**4 * emitx * emity * sigma_e * bunch_length
        )
        integrals = _bjorken_mtingwa(optics, gamma, emitx, emity, sigma_e)
    if model == "Nagaitsev":
        integrals /= np.array([emitx, emity, sigma_e**2])
    return cst * integrals


def _npart(ring: Lattice, bunch_current):
    """Number of particles in bunches of given current."""
    return bunch_current / (ring.revolution_frequency * qe * abs(ring.particle.charge))


def ibs_rates(ring: Lattice, emitx: float, emity: float, sigma_e: float,
              bunch_length: float, bunch_current: float | None = None, *,
              model: str = "Nagaitsev", step: float = 0.1) -> np.ndarray:
    r"""IBS emittance growth rates.

    Parameters:
        ring:           Lattice description
        emitx:          Horizontal emittance [m]
        emity:          Vertical emittance [m]
        sigma_e:        Relative momentum spread
        bunch_length:   RMS bunch length [m]
        bunch_current:  Bunch current [A]. Default: current of the 1st bunch
          of *ring*
        model:          ``"Nagaitsev"`` or ``"Bjorken-Mtingwa"`` (also
          ``"B&M"``), case-insensitive
        step:           Maximum distance between optics sampling points [m]

    Returns:
        rates:  Emittance growth rates :math:`\frac{1}{\epsilon}
          \frac{d\epsilon}{dt}` [1/s] in the horizontal, vertical and
          longitudinal planes. The longitudinal rate is the one of
          :math:`\sigma_\delta^2`. Amplitude growth rates, as returned by
          xsuite, are half of these values.
    """
    model = _check_model(model)
    if bunch_current is None:
        bunch_current = ring.bunch_currents[0]
    return _growth_rates(_ibs_optics(ring, step), model, _beam_constants(ring),
                         _npart(ring, bunch_current), emitx, emity, sigma_e,
                         bunch_length)


def ibs_equilibrium(ring: Lattice, bunch_current: float | None = None, *,
                    model: str = "Nagaitsev", coupling: float = 0.0,
                    constraint: str = "coupling", rtol: float = 1.0e-6,
                    max_steps: int = 100000, step: float = 0.1) -> dict:
    r"""Equilibrium emittances with IBS and synchrotron radiation.

    Integrates :math:`\frac{d\epsilon}{dt} = -\frac{2}{\tau}(\epsilon -
    \epsilon_0) + r_{IBS}\,\epsilon` in the 3 planes, starting from the
    radiation equilibrium :math:`\epsilon_0` and with adaptive time steps,
    until the relative change of the emittances is below *rtol*. The ratio
    of bunch length to momentum spread is kept constant.

    Parameters:
        ring:           Lattice description. Radiation parameters are
          computed from the radiation integrals
        bunch_current:  Bunch current [A]. Default: current of the 1st bunch
          of *ring*
        model:          ``"Nagaitsev"`` or ``"Bjorken-Mtingwa"`` (also
          ``"B&M"``), case-insensitive
        coupling:       Emittance ratio :math:`\kappa=\epsilon_y/\epsilon_x`.
          Default: 0, use the vertical equilibrium emittance of *ring*
        constraint:     How the vertical emittance follows the horizontal
          one when *coupling* is non-zero:

          * ``"coupling"``: betatron coupling, :math:`\epsilon_x +
            \epsilon_y` is shared according to :math:`\kappa`,
          * ``"excitation"``: vertical excitation, :math:`\epsilon_y =
            \kappa\,\epsilon_x`.
        rtol:           Relative tolerance for convergence
        max_steps:      Maximum number of time steps
        step:           Maximum distance between optics sampling points [m]

    Returns:
        result: dictionary with keys:

          * ``emittances``: equilibrium emittances [m], horizontal, vertical
            and longitudinal (:math:`\sigma_\delta\sigma_z`),
          * ``sigma_e``: equilibrium momentum spread,
          * ``bunch_length``: equilibrium bunch length [m],
          * ``rates``: IBS emittance growth rates at equilibrium [1/s],
          * ``time``: (nsteps,) time [s],
          * ``history``: (nsteps, 3) evolution of the emittances [m].
    """
    model = _check_model(model)
    constraint = constraint.lower()
    if constraint not in ("coupling", "excitation"):
        raise ValueError("constraint must be 'coupling' or 'excitation'")
    if bunch_current is None:
        bunch_current = ring.bunch_currents[0]
    optics = _ibs_optics(ring, step)
    consts = _beam_constants(ring)
    npart = _npart(ring, bunch_current)

    rp = ring.disable_6d(copy=True).radiation_parameters()
    damping = 2.0 / rp.Tau  # emittance damping rates
    ratio = rp.sigma_l / rp.sigma_e
    eq0 = np.nan_to_num(np.array([*rp.emittances[:2], rp.sigma_l * rp.sigma_e]))
    if coupling > 0.0 and constraint == "coupling":
        jx, jy = rp.J[:2]
        eq0[0] /= 1.0 + coupling * jy / jx
        eq0[1] = coupling * eq0[0]
    elif coupling > 0.0:
        eq0[1] = coupling * eq0[0]
    if eq0[1] <= 0.0:
        raise ValueError("No vertical emittance: set the coupling")

    def constrain(emit):
        if coupling > 0.0 and constraint == "coupling":
            emit[0] = (emit[0] + emit[1]) / (1.0 + coupling)
            emit[1] = coupling * emit[0]
        elif coupling > 0.0:
            emit[1] = coupling * emit[0]
        return emit

    emit = eq0.copy()
    time, history = [0.0], [emit.copy()]
    for _ in range(max_steps):
        rates = _growth_rates(optics, model, consts, npart, emit[0], emit[1],
                              np.sqrt(emit[2] / ratio), np.sqrt(emit[2] * ratio))
        # Time step as in xsuite: 1% of the fastest amplitude rate
        dt = 0.02 / np.max(np.abs(np.concatenate((rates, damping))))
        new = constrain(emit + dt * (-damping * (emit - eq0) + rates * emit))
        time.append(time[-1] + dt)
        history.append(new.copy())
        converged = np.max(np.abs(new / emit - 1.0)) < rtol
        emit = new
        if converged:
            break
    else:
        raise RuntimeError(f"IBS equilibrium not reached in {max_steps} steps")

    return {
        "emittances": emit,
        "sigma_e": np.sqrt(emit[2] / ratio),
        "bunch_length": np.sqrt(emit[2] * ratio),
        "rates": rates,
        "time": np.array(time),
        "history": np.array(history),
    }


class IBSElement(Collective, Element):
    """Element applying IBS momentum kicks during tracking.

    Every turn, each particle receives random kicks in :math:`p_x`,
    :math:`p_y` and :math:`\\delta`. The kick amplitudes are computed from
    analytical growth rates (see :py:func:`ibs_rates`) evaluated with the
    emittances, momentum spread and bunch length of the tracked bunch. They
    are weighted by the local longitudinal line density.

    The momentum kick is compensated in :math:`x, p_x, y, p_y` by the local
    dispersion, so that it does not change the betatron coordinates: the
    dispersive contribution to transverse growth is already included in the
    growth rates.

    Multi-bunch beams are handled according to the filling pattern of *ring*
    at the time the element is created.

    Random kicks use :py:obj:`at.random.thread <.random>`, seeded with
    :py:meth:`at.random.reset <.random.reset>`.
    """

    default_pass = {False: "IdentityPass", True: "pyIBSPass"}

    def __init__(self, family_name: str, ring: Lattice, *,
                 model: str = "Nagaitsev", update_turns: int = 100,
                 nslice: int = 51, step: float = 0.1, refpt: int = 0,
                 **kwargs):
        """
        Parameters:
            family_name:    Element name
            ring:           Lattice in which the element will be inserted
            model:          Growth rate model: ``"Nagaitsev"`` or
              ``"Bjorken-Mtingwa"`` (also ``"B&M"``), case-insensitive
            update_turns:   Number of turns between updates of the growth
              rates
            nslice:         Number of bins of the longitudinal line density.
              Use 1 for uniform kicks
            step:           Maximum distance between optics sampling points
              [m]
            refpt:          Index in *ring* where the element will be
              inserted. Default: 0, entrance (or end) of *ring*.

        Attributes:
            rates:  Emittance growth rates [1/s] of each bunch, shape
              (nbunch, 3), computed at the last update
        """
        kwargs.setdefault("PassMethod", self.default_pass[True])
        self.model = _check_model(model)
        self.update_turns = int(update_turns)
        self.nslice = int(nslice)
        self._consts = _beam_constants(ring)
        # Time between passes: the lattice may be one cell of the ring
        self._dt = 1.0 / ring.cell_revolution_frequency
        self._npart = _npart(ring, ring.bunch_currents)
        self._optics = _ibs_optics(ring, step)
        _, _, ld = ring.disable_6d(copy=True).get_optics(refpts=refpt)
        self._local_beta = ld.beta[0]
        self._local_disp = ld.dispersion[0]
        self.clear_history()
        super().__init__(family_name, **kwargs)

    def clear_history(self):
        """Reset the turn counter and growth rates."""
        self._turn = 0
        self.rates = np.zeros((len(self._npart), 3))

    def __repr__(self):
        """Simplified __repr__: the element cannot be rebuilt without *ring*"""
        att = {k: v for (k, v) in self.items() if not k.startswith("_")}
        return f"{self.__class__.__name__}({att})"

    def _beam_parameters(self, bunch):
        """Emittances, momentum spread and bunch length of a bunch."""
        delta = bunch[4]
        betatron = bunch[:4] - np.outer(self._local_disp, delta)
        emitx = np.sqrt(np.linalg.det(np.cov(betatron[:2])))
        emity = np.sqrt(np.linalg.det(np.cov(betatron[2:])))
        return emitx, emity, np.std(delta), np.std(bunch[5])

    def _line_density(self, ct):
        """Line density at each particle position, normalised to mean 1."""
        if self.nslice <= 1:
            return np.ones_like(ct)
        hist, edges = np.histogram(ct, bins=self.nslice)
        rho = hist[np.clip(np.digitize(ct, edges) - 1, 0, self.nslice - 1)]
        return rho / rho.mean()

    def kick(self, bunch: np.ndarray, ibunch: int = 0) -> None:
        """Apply the IBS kicks to the alive particles of a bunch, in place.

        Parameters:
            bunch:  (6, N) particle coordinates
            ibunch: Bunch index in the filling pattern
        """
        alive = np.flatnonzero(~np.isnan(bunch[0]))
        if len(alive) < 2:
            return
        part = bunch[:, alive]
        emitx, emity, sigma_e, sigma_s = self._beam_parameters(part)
        if self._turn % self.update_turns == 0:
            self.rates[ibunch] = _growth_rates(
                self._optics, self.model, self._consts,
                self._npart[ibunch], emitx, emity, sigma_e, sigma_s,
            )
        # <dp^2> per pass = 2 * d(emittance) / beta for the transverse planes
        rx, ry, rp = np.maximum(self.rates[ibunch], 0.0) * self._dt
        bx, by = self._local_beta
        sigma = np.array([np.sqrt(2.0 * rx * emitx / bx),
                          np.sqrt(2.0 * ry * emity / by),
                          np.sqrt(2.0 * rp) * sigma_e])
        rnd = random.thread.standard_normal((3, len(alive)))
        dp = sigma[:, None] * rnd * np.sqrt(self._line_density(part[5]))
        part[[1, 3, 4]] += dp
        part[:4] += np.outer(self._local_disp, dp[2])
        bunch[:, alive] = part

    def track_turn(self, rin: np.ndarray) -> None:
        """Apply the IBS kicks to all bunches and increment the turn counter."""
        nbunch = len(self._npart)
        for ib in range(nbunch):
            bunch = rin[:, ib::nbunch]
            self.kick(bunch, ib)
        self._turn += 1
