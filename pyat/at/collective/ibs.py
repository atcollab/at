"""Intra-beam scattering (IBS).

Analytical IBS growth rates, equilibrium emittances with radiation damping,
and a collective element applying the corresponding random momentum kicks
during tracking.

The growth rates follow the Bjorken-Mtingwa model, J.D. Bjorken,
S.K. Mtingwa, Part. Accel. 13, 115 (1983), as modified in MAD-X to include
the vertical dispersion: F. Antoniou, F. Zimmermann, CERN-ATS-2012-066. The
Coulomb logarithm is computed as in MAD-X (``twclog``). They are computed in
C, by the same code as the ``IBSPass`` pass method.

The local momentum distribution entering the Bjorken-Mtingwa integrals is
built from the 6D mode matrices of :py:func:`.linopt6`, so that betatron
and synchro-betatron coupling are taken into account, as in K. Kubo,
K. Oide, PRSTAB 4, 124401 (2001). The emittances are the mode emittances.
For an uncoupled lattice this reduces to the MAD-X formulas.

The kick model follows R. Bruce et al., PRSTAB 13, 091001 (2010), as also
implemented in mbtrack2, xfields and elegant.
"""

from __future__ import annotations

__all__ = ["IBSConstraint", "IBSElement", "ibs_equilibrium", "ibs_rates"]

from enum import Enum

import numpy as np
from scipy.optimize import root

from ..constants import e_mass, qe
from ..lattice import All, AtError, Collective, Element, Lattice, check_6d

# noinspection PyProtectedMember
from ..lattice.elements.conversions import _array
from ..physics.amat import jmat
from ._ibs import growth_rates


class IBSConstraint(Enum):
    r"""Evolution of the vertical emittance in :py:func:`ibs_equilibrium`."""

    #: Independent vertical emittance, with its own damping and IBS rate
    UNCOUPLED = 1
    #: Betatron coupling: :math:`\epsilon_x + \epsilon_y` is shared with a
    #: constant ratio :math:`\kappa=\epsilon_y/\epsilon_x`. The radiation
    #: equilibrium :math:`\epsilon_{x0}` is reduced so that
    #: :math:`\epsilon_{x0} + \frac{J_y}{J_x}\,\epsilon_{y0}` is the natural
    #: emittance, and :math:`\kappa=\epsilon_{y0}/\epsilon_{x0}`
    COUPLED = 2
    #: Vertical excitation: :math:`\epsilon_y = \kappa\,\epsilon_x` with
    #: :math:`\kappa=\epsilon_{y0}/\epsilon_{x0}`
    FIXED_RATIO = 3


_J = jmat(3)
_MOM = [1, 4, 3]


@check_6d(True)
def _ibs_optics(ring: Lattice, step: float) -> np.ndarray:
    r"""Optics sampled around the ring.

    Returns:
        optics: (21, N) array: trapezoidal integration weights normalised to
          1, for the ring averages of the growth rates and Coulomb logarithm,
          :math:`D_x, D_y`, and the upper triangles of the momentum blocks
          :math:`P_k` of :math:`J^TR_kJ` in the Bjorken-Mtingwa coordinates
          :math:`(p_x, \delta/\gamma, p_y)`, for the 3 modes. The longitudinal
          mode matrix is divided by its :math:`\delta\delta` element at the
          entrance of the ring, so that its emittance is
          :math:`\sigma_\delta^2`
    """
    _, _, ld = ring.slice(size=step).get_optics(refpts=All)
    w = np.convolve(np.diff(ld.s_pos), [1.0, 1.0])
    q = _J.T @ ld.R @ _J
    q[:, 2] *= ld.R[0, 2, 4, 4]
    scale = np.array([1.0, ring.gamma, 1.0])
    pmat = q[:, :, _MOM][:, :, :, _MOM] * np.outer(scale, scale)
    iu = np.triu_indices(3)
    return np.vstack(
        (
            w / w.sum(),
            ld.dispersion[:, [0, 2]].T,
            pmat[:, :, iu[0], iu[1]].reshape(len(w), 18).T,
        )
    )


def _beam_constants(ring: Lattice) -> np.ndarray:
    """Energy [eV], rest energy [eV], charge and relativistic beta.

    As in :py:class:`.Lattice`, a particle with zero rest energy
    (``"relativistic"``) has the electron mass and beta = 1."""
    particle = ring.particle
    mass = particle.rest_energy if particle.rest_energy > 0.0 else e_mass
    return np.array([ring.energy, mass, particle.charge, ring.beta])


def _npart(ring: Lattice, bunch_current):
    """Number of particles in bunches of given current."""
    return bunch_current / (ring.revolution_frequency * qe * abs(ring.particle.charge))


@check_6d(True)
def ibs_rates(
    ring: Lattice,
    emitx: float,
    emity: float,
    sigma_e: float,
    bunch_length: float,
    bunch_current: float | None = None,
    *,
    step: float = 0.1,
) -> np.ndarray:
    r"""IBS emittance growth rates.

    Parameters:
        ring:           Lattice description, must be 6D
        emitx:          Horizontal (mode 1) emittance [m]
        emity:          Vertical (mode 2) emittance [m]
        sigma_e:        Relative momentum spread
        bunch_length:   RMS bunch length [m]
        bunch_current:  Bunch current [A]. Default: current of the 1st bunch
          of *ring*
        step:           Maximum distance between optics sampling points [m]

    Returns:
        rates:  Emittance growth rates :math:`\frac{1}{\epsilon}
          \frac{d\epsilon}{dt}` [1/s] of the horizontal, vertical and
          longitudinal modes. The longitudinal emittance is scaled to
          :math:`\sigma_\delta^2`. Amplitude growth rates, as returned by
          xsuite, are half of these values.
    """
    if bunch_current is None:
        bunch_current = ring.bunch_currents[0]
    return growth_rates(
        _ibs_optics(ring, step),
        *_beam_constants(ring),
        _npart(ring, bunch_current),
        emitx,
        emity,
        sigma_e,
        bunch_length,
    )


@check_6d(True)
def ibs_equilibrium(
    ring: Lattice,
    bunch_current: float | None = None,
    *,
    emity: float | None = None,
    constraint: IBSConstraint = IBSConstraint.UNCOUPLED,
    rtol: float = 1.0e-6,
    max_steps: int = 100000,
    step: float = 0.1,
    history: bool = False,
) -> dict:
    r"""Equilibrium emittances with IBS and synchrotron radiation.

    Solves for the steady state of :math:`\frac{d\epsilon}{dt} =
    -\frac{2}{\tau}(\epsilon - \epsilon_0) + r_{IBS}\,\epsilon` in the
    3 planes, where :math:`\epsilon` are the mode emittances and
    :math:`\epsilon_0` the radiation equilibrium. The momentum spread and
    bunch length scale as :math:`\sqrt{\epsilon_z/\epsilon_{z0}}`.

    The steady state is found with a root solver. If it fails, or if
    *history* is :py:obj:`True`, the equations are instead integrated in
    time from :math:`\epsilon_0` with adaptive time steps, until
    :math:`|\frac{d\epsilon}{dt}| / (\frac{2}{\tau}\epsilon)` is below *rtol*
    in all planes, the same criterion as for the root solver.

    Parameters:
        ring:           Lattice description, must be 6D. Radiation
          parameters are computed with :py:func:`.envelope_parameters`
        bunch_current:  Bunch current [A]. Default: current of the 1st bunch
          of *ring*
        emity:          Vertical equilibrium emittance without IBS [m].
          Default: vertical equilibrium emittance of *ring*
        constraint:     How the vertical emittance evolves with IBS, see
          :py:class:`IBSConstraint`. Default:
          :py:attr:`IBSConstraint.UNCOUPLED`
        rtol:           Relative tolerance for convergence
        max_steps:      Maximum number of time steps
        step:           Maximum distance between optics sampling points [m]
        history:        Integrate in time and return the evolution of the
          emittances

    Returns:
        result: dictionary with keys:

          * ``emittances``: equilibrium mode emittances [m], horizontal,
            vertical and longitudinal,
          * ``sigma_e``: equilibrium momentum spread,
          * ``bunch_length``: equilibrium bunch length [m],
          * ``rates``: IBS emittance growth rates at equilibrium [1/s],
          * ``time``: (nsteps,) time [s], only with time integration,
          * ``history``: (nsteps, 3) evolution of the emittances [m], only
            with time integration.
    """
    if bunch_current is None:
        bunch_current = ring.bunch_currents[0]
    optics = _ibs_optics(ring, step)
    consts = _beam_constants(ring)
    npart = _npart(ring, bunch_current)

    rp = ring.envelope_parameters()
    damping = 2.0 / rp.Tau
    eq0 = rp.emittances.copy()
    if emity is not None:
        eq0[1] = emity
    if not eq0[1] > 1.0e-14:
        msg = "No vertical emittance: set emity"
        raise AtError(msg)
    if constraint is IBSConstraint.COUPLED:
        jx, jy = rp.J[:2]
        enat = rp.emittances[0] + jy / jx * rp.emittances[1]
        emax = enat * jx / (jx + jy)
        if eq0[1] > emax:
            msg = (
                f"emity larger than emitx: the maximum, at full coupling, "
                f"is {emax:.4g} m"
            )
            raise AtError(msg)
        eq0[0] = enat - jy / jx * eq0[1]
    kappa = eq0[1] / eq0[0]
    tied = constraint is not IBSConstraint.UNCOUPLED
    free = [0, 2] if tied else [0, 1, 2]

    def long_sizes(emit):
        scale = np.sqrt(emit[2] / eq0[2])
        return rp.sigma_e * scale, rp.sigma_l * scale

    def constrain(emit):
        if constraint is IBSConstraint.COUPLED:
            emit[0] = (emit[0] + emit[1]) / (1.0 + kappa)
            emit[1] = kappa * emit[0]
        elif constraint is IBSConstraint.FIXED_RATIO:
            emit[1] = kappa * emit[0]
        return emit

    def ibs(emit):
        return growth_rates(optics, *consts, npart, emit[0], emit[1], *long_sizes(emit))

    def derivative(emit, rates):
        return -damping * (emit - eq0) + rates * emit

    def normalized(deriv, emit):
        res = deriv / (damping * emit)
        if constraint is IBSConstraint.COUPLED:
            res[0] += deriv[1] / (damping[0] * emit[0])
        return res[free]

    def result(emit, **kwargs):
        sigma_e, bunch_length = long_sizes(emit)
        return dict(
            emittances=emit,
            sigma_e=sigma_e,
            bunch_length=bunch_length,
            rates=ibs(emit),
            **kwargs,
        )

    def emittances(u):
        emit = eq0.copy()
        emit[free] = np.exp(u)
        if tied:
            emit[1] = kappa * emit[0]
        return emit

    def residual(emit):
        return normalized(derivative(emit, ibs(emit)), emit)

    if not history:
        sol = root(lambda u: residual(emittances(u)), np.log(eq0[free]), tol=rtol)
        if sol.success and np.max(np.abs(sol.fun)) < rtol:
            return result(emittances(sol.x))

    emit = eq0
    times, evolution = [0.0], [emit]
    for _ in range(max_steps):
        rates = ibs(emit)
        deriv = derivative(emit, rates)
        if np.max(np.abs(normalized(deriv, emit))) < rtol:
            break
        dt = 0.02 / np.max(np.abs(np.concatenate((rates, damping))))
        emit = constrain(emit + dt * deriv)
        times.append(times[-1] + dt)
        evolution.append(emit)
    else:
        msg = f"IBS equilibrium not reached in {max_steps} steps"
        raise AtError(msg)

    return result(emit, time=np.array(times), history=np.array(evolution))


class IBSElement(Collective, Element):
    """Element applying IBS momentum kicks during tracking.

    Every turn, each particle receives a random 6D kick. The kick amplitudes
    are computed from analytical growth rates (see :py:func:`ibs_rates`)
    evaluated with the transverse mode emittances, momentum spread and bunch
    length of the tracked bunch. They are weighted by the local longitudinal
    line density.

    The transverse emittances are the projections
    :math:`\\frac{1}{2}\\mathrm{Tr}(J^TR_kJ\\Sigma)` of the sigma matrix of
    the bunch on the modes :math:`R_k` of the lattice at the element: for a
    matched bunch they are the eigen emittances, for a mismatched one the
    emittances it filaments to. The element must therefore be created with
    the lattice in which it is tracked. IBS is not defined for a bunch with
    zero emittance or spread in a plane: such a bunch receives no kick and
    its growth rates are set to NaN. The kick of each mode has the shape
    :math:`R_k` of the mode at the element, so that each mode grows at its
    own rate without exciting the others or mismatching the bunch.

    The kicks are computed by the ``IBSPass`` C pass method, with random
    numbers from the AT C generators, seeded with :py:func:`.reset_rng`.
    The number of particles of each bunch is computed from the bunch
    currents of the tracked lattice.
    """

    default_pass = {False: "IdentityPass", True: "IBSPass"}
    _conversions = dict(
        Element._conversions,
        UpdateTurns=int,
        NSlice=int,
        _npart_per_amp=float,
        _optics=lambda v: _array(v, (21, -1)),
        _local_kicks=lambda v: _array(v, (36,)),
        _local_modes=lambda v: _array(v, (72,)),
        _rates=lambda v: _array(v, (-1, 3)),
    )

    def __init__(
        self,
        family_name: str,
        ring: Lattice,
        *,
        update_turns: int = 1,
        nslice: int = 51,
        step: float = 0.1,
        refpt: int = 0,
        **kwargs,
    ):
        """
        Parameters:
            family_name:    Element name
            ring:           Lattice in which the element will be inserted,
              must be 6D
            update_turns:   Number of turns between updates of the growth
              rates. Between updates, the kicks are rescaled with the
              current beam parameters at constant rates: the interval must
              be much shorter than the IBS growth time, otherwise the
              emittances grow exponentially. Default: 1, the update costs
              little compared with the tracking
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
        self._npart_per_amp = _npart(ring, 1.0)
        self._optics = _ibs_optics(ring, step)
        _, _, ld = ring.get_optics(refpts=refpt)
        scale = np.sqrt([1.0, 1.0, ld.R[0, 2, 4, 4]])
        a = ld.A[0].reshape(6, 3, 2).transpose(1, 0, 2) / scale[:, None, None]
        self._local_kicks = a.ravel()
        self._local_modes = (_J.T @ ld.R[0, :2] @ _J).ravel()
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
        """Simplified __repr__: the element cannot be rebuilt without *ring*."""
        att = {k: v for (k, v) in self.items() if not k.startswith("_")}
        return f"{self.__class__.__name__}({att})"
