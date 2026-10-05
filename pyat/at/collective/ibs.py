"""Intra-beam scattering (IBS).

Analytical IBS growth rates, equilibrium emittances with radiation damping,
and a collective element applying the corresponding random momentum kicks
during tracking.

The growth rates follow the Bjorken-Mtingwa model, J.D. Bjorken,
S.K. Mtingwa, Part. Accel. 13, 115 (1983), as modified in MAD-X to include
the vertical dispersion: F. Antoniou, F. Zimmermann, CERN-ATS-2012-066. The
Coulomb logarithm is computed as in MAD-X (``twclog``). They are computed in
C, by the same code as the ``IBSPass`` pass method.

The kick model follows R. Bruce et al., PRSTAB 13, 091001 (2010), as also
implemented in mbtrack2, xfields and elegant.
"""

from __future__ import annotations

__all__ = ["ibs_rates", "ibs_equilibrium", "IBSElement"]

import numpy as np
from scipy.optimize import root

from ..constants import e_mass, qe
from ..lattice import All, Collective, Element, Lattice
# noinspection PyProtectedMember
from ..lattice.elements.conversions import _array
from ._ibs import growth_rates

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
    """Energy [eV], rest energy [eV], charge and relativistic beta.

    As in :py:class:`.Lattice`, a particle with zero rest energy
    (``"relativistic"``) has the electron mass and beta = 1."""
    particle = ring.particle
    mass = particle.rest_energy if particle.rest_energy > 0.0 else e_mass
    return np.array([ring.energy, mass, particle.charge, ring.beta])


def _growth_rates(optics: np.ndarray, consts: np.ndarray, npart: float,
                  emitx: float, emity: float, sigma_e: float,
                  bunch_length: float) -> np.ndarray:
    """Emittance growth rates [1/s] for given beam parameters, computed in C
    as in the IBSPass pass method."""
    return growth_rates(optics, *consts, npart, emitx, emity, sigma_e,
                        bunch_length)


def _npart(ring: Lattice, bunch_current):
    """Number of particles in bunches of given current."""
    return bunch_current / (ring.revolution_frequency * qe * abs(ring.particle.charge))


def ibs_rates(ring: Lattice, emitx: float, emity: float, sigma_e: float,
              bunch_length: float, bunch_current: float | None = None, *,
              step: float = 0.1) -> np.ndarray:
    r"""IBS emittance growth rates.

    Parameters:
        ring:           Lattice description
        emitx:          Horizontal emittance [m]
        emity:          Vertical emittance [m]
        sigma_e:        Relative momentum spread
        bunch_length:   RMS bunch length [m]
        bunch_current:  Bunch current [A]. Default: current of the 1st bunch
          of *ring*
        step:           Maximum distance between optics sampling points [m]

    Returns:
        rates:  Emittance growth rates :math:`\frac{1}{\epsilon}
          \frac{d\epsilon}{dt}` [1/s] in the horizontal, vertical and
          longitudinal planes. The longitudinal rate is the one of
          :math:`\sigma_\delta^2`. Amplitude growth rates, as returned by
          xsuite, are half of these values.
    """
    if bunch_current is None:
        bunch_current = ring.bunch_currents[0]
    return _growth_rates(_ibs_optics(ring, step), _beam_constants(ring),
                         _npart(ring, bunch_current), emitx, emity, sigma_e,
                         bunch_length)


def ibs_equilibrium(ring: Lattice, bunch_current: float | None = None, *,
                    coupling: float = 0.0,
                    constraint: str = "coupling", rtol: float = 1.0e-6,
                    max_steps: int = 100000, step: float = 0.1,
                    history: bool = False) -> dict:
    r"""Equilibrium emittances with IBS and synchrotron radiation.

    Solves for the steady state of :math:`\frac{d\epsilon}{dt} =
    -\frac{2}{\tau}(\epsilon - \epsilon_0) + r_{IBS}\,\epsilon` in the
    3 planes, where :math:`\epsilon_0` is the radiation equilibrium. The
    ratio of bunch length to momentum spread is kept constant.

    The steady state is found with a root solver. If it fails, or if
    *history* is :py:obj:`True`, the equations are instead integrated in
    time from :math:`\epsilon_0` with adaptive time steps, until the
    relative change of the emittances is below *rtol*.

    Parameters:
        ring:           Lattice description. Radiation parameters are
          computed from the radiation integrals
        bunch_current:  Bunch current [A]. Default: current of the 1st bunch
          of *ring*
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
        history:        Integrate in time and return the evolution of the
          emittances

    Returns:
        result: dictionary with keys:

          * ``emittances``: equilibrium emittances [m], horizontal, vertical
            and longitudinal (:math:`\sigma_\delta\sigma_z`),
          * ``sigma_e``: equilibrium momentum spread,
          * ``bunch_length``: equilibrium bunch length [m],
          * ``rates``: IBS emittance growth rates at equilibrium [1/s],
          * ``time``: (nsteps,) time [s], only with time integration,
          * ``history``: (nsteps, 3) evolution of the emittances [m], only
            with time integration.
    """
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

    def ibs(emit):
        return _growth_rates(optics, consts, npart, emit[0], emit[1],
                             np.sqrt(emit[2] / ratio), np.sqrt(emit[2] * ratio))

    def derivative(emit):
        return -damping * (emit - eq0) + ibs(emit) * emit

    def result(emit, **kwargs):
        return dict(emittances=emit, sigma_e=np.sqrt(emit[2] / ratio),
                    bunch_length=np.sqrt(emit[2] * ratio), rates=ibs(emit),
                    **kwargs)

    # Unknowns: log of the independent emittances. With coupling, the
    # vertical emittance follows the horizontal one
    def emittances(u):
        emit = np.exp(u)
        if coupling > 0.0:
            emit = np.array([emit[0], coupling * emit[0], emit[1]])
        return emit

    def residual(u):
        emit = emittances(u)
        deriv = derivative(emit)
        if coupling > 0.0 and constraint == "coupling":
            deriv[0] += deriv[1]
        deriv /= damping * emit
        return deriv[[0, 2]] if coupling > 0.0 else deriv

    if not history:
        u0 = np.log(eq0[[0, 2]] if coupling > 0.0 else eq0)
        sol = root(residual, u0, tol=rtol)
        if sol.success and np.max(np.abs(residual(sol.x))) < rtol:
            return result(emittances(sol.x))

    emit = eq0.copy()
    times, evolution = [0.0], [emit.copy()]
    for _ in range(max_steps):
        rates = ibs(emit)
        # Time step as in xsuite: 1% of the fastest amplitude rate
        dt = 0.02 / np.max(np.abs(np.concatenate((rates, damping))))
        new = constrain(emit + dt * (-damping * (emit - eq0) + rates * emit))
        times.append(times[-1] + dt)
        evolution.append(new.copy())
        converged = np.max(np.abs(new / emit - 1.0)) < rtol
        emit = new
        if converged:
            break
    else:
        raise RuntimeError(f"IBS equilibrium not reached in {max_steps} steps")

    return result(emit, time=np.array(times), history=np.array(evolution))


class IBSElement(Collective, Element):
    """Element applying IBS momentum kicks during tracking.

    Every turn, each particle receives random kicks in :math:`p_x`,
    :math:`p_y` and :math:`\\delta`. The kick amplitudes are computed from
    analytical growth rates (see :py:func:`ibs_rates`) evaluated with the
    emittances, momentum spread and bunch length of the tracked bunch. They
    are weighted by the local longitudinal line density.

    The emittances are the eigen emittances of the 6D sigma matrix of the
    bunch, as in :py:func:`.emittances_from_beam`. IBS is not defined for a
    bunch with a singular sigma matrix (for instance a plane with zero
    spread): such a bunch receives no kick and its growth rates are set to
    NaN. The momentum kick is compensated in
    :math:`x, p_x, y, p_y` by the local dispersion, so that it does not
    change the betatron coordinates: the dispersive contribution to
    transverse growth is already included in the growth rates.

    The kicks are computed by the ``IBSPass`` C pass method, with random
    numbers from the AT C generators, seeded with :py:func:`.reset_rng`.
    The number of particles of each bunch is computed from the bunch
    currents of the tracked lattice.
    """

    default_pass = {False: "IdentityPass", True: "IBSPass"}
    _conversions = dict(Element._conversions, UpdateTurns=int,
                        NSlice=int, _trev=float,
                        _optics=lambda v: _array(v, (9, -1)),
                        _local_beta=lambda v: _array(v, (2,)),
                        _local_disp=lambda v: _array(v, (4,)),
                        _rates=lambda v: _array(v, (-1, 3)))

    def __init__(self, family_name: str, ring: Lattice, *,
                 update_turns: int = 100, nslice: int = 51,
                 step: float = 0.1, refpt: int = 0,
                 **kwargs):
        """
        Parameters:
            family_name:    Element name
            ring:           Lattice in which the element will be inserted
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
        self.UpdateTurns = update_turns
        self.NSlice = nslice
        self._trev = 1.0 / ring.revolution_frequency
        self._optics = _ibs_optics(ring, step)
        _, _, ld = ring.disable_6d(copy=True).get_optics(refpts=refpt)
        self._local_beta = ld.beta[0]
        self._local_disp = ld.dispersion[0]
        self._rates = np.full((ring.nbunch, 3), np.nan)
        super().__init__(family_name, **kwargs)

    @property
    def rates(self) -> np.ndarray:
        """Emittance growth rates [1/s] of each bunch, shape (nbunch, 3)."""
        return self._rates

    def clear_history(self, ring: Lattice | None = None):
        """Reset the growth rates: they are recomputed at the next pass.

        Parameters:
            ring:   If given, the rates are resized for the bunches of *ring*
        """
        nbunch = len(self._rates) if ring is None else ring.nbunch
        self._rates = np.full((nbunch, 3), np.nan)

    def __repr__(self):
        """Simplified __repr__: the element cannot be rebuilt without *ring*"""
        att = {k: v for (k, v) in self.items() if not k.startswith("_")}
        return f"{self.__class__.__name__}({att})"
