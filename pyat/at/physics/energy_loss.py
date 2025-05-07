from __future__ import annotations

__all__ = ["get_energy_loss", "set_cavity_phase", "ELossMethod", "get_timelag_fromU0"]

from enum import Enum
from warnings import warn
from collections.abc import Sequence

import numpy as np
from scipy.optimize import least_squares

from at.constants import clight, Cgamma
from at.lattice import Lattice, Dipole, Wiggler, RFCavity, Refpts, EnergyLoss
from at.lattice import check_radiation, AtError, AtWarning
from at.lattice import get_bool_index, set_value_refpts
from at.lattice import DConstant


class ELossMethod(Enum):
    """methods for the computation of energy losses"""

    #: The losses are obtained from
    #: :math:`E_{loss}=C_\gamma/2\pi . E^4 . I_2`.
    #: Takes into account bending magnets and wigglers.
    INTEGRAL = 1
    #: The losses are obtained by tracking without cavities.
    #: Needs radiation ON, takes into account all radiating elements
    TRACKING = 2


def get_energy_loss(
    ring: Lattice, method: ELossMethod | None = ELossMethod.INTEGRAL
) -> float:
    """Computes the energy loss per turn

    Parameters:
        ring:           Lattice description
        method:         Method for energy loss computation.
          See :py:class:`ELossMethod`.

    Returns:
        eloss (float):  Energy loss per turn [eV]
    """

    # noinspection PyShadowingNames
    def integral(ring):
        """Losses = Cgamma / 2pi * EGeV^4 * i2"""

        def wiggler_i2(wiggler: Wiggler):
            rhoinv = wiggler.Bmax / ring.BRho
            coefh = wiggler.By[1, :]
            coefv = wiggler.Bx[1, :]
            return (
                wiggler.Length
                * (np.sum(coefh * coefh) + np.sum(coefv * coefv))
                * rhoinv**2
                / 2
            )

        def dipole_i2(dipole: Dipole):
            return dipole.BendingAngle**2 / dipole.Length

        def eloss_i2(eloss: EnergyLoss):
            return eloss.EnergyLoss / coef

        i2 = 0.0
        coef = Cgamma / 2.0 / np.pi * ring.energy**4
        for el in ring:
            if isinstance(el, Dipole):
                i2 += dipole_i2(el)
            elif isinstance(el, Wiggler) and el.PassMethod != "DriftPass":
                i2 += wiggler_i2(el)
            elif isinstance(el, EnergyLoss) and el.PassMethod != "IdentityPass":
                i2 += eloss_i2(el)
        e_loss = coef * i2
        return e_loss

    # noinspection PyShadowingNames
    @check_radiation(True)
    def tracking(ring):
        """Losses from tracking"""
        energy = ring.energy
        particle = ring.particle
        delta = 0.0
        for e in ring:
            if e.PassMethod.endswith("RadPass"):
                ot = e.track(np.zeros(6), energy=energy, particle=particle)
                delta += ot[4]
        return -delta * energy

    if isinstance(method, str):
        method = ELossMethod[method.upper()]
        warn(FutureWarning(f"You should use {method!s}"), stacklevel=2)
    if method is ELossMethod.INTEGRAL:
        return ring.periodicity * integral(ring)
    elif method == ELossMethod.TRACKING:
        return ring.periodicity * tracking(ring)
    else:
        raise AtError(f"Invalid method: {method}")


# noinspection PyPep8Naming
def get_timelag_fromU0(
    ring: Lattice,
    *,
    method: ELossMethod | None = ELossMethod.TRACKING,
    cavpts: Refpts | None = None,
    divider: int | None = 4,
    ts_tol: float | None = None,
) -> tuple[Sequence[float], float]:
    """
    Get the TimeLag attribute of RF cavities based on frequency,
    voltage and energy loss per turn, so that the synchronous phase is zero.
    Used in set_cavity_phase()

    Parameters:
        ring:               Lattice description
        method:             Method for energy loss computation.
          See :py:class:`ELossMethod`.
        cavpts:             Cavity location. If None, use all cavities.
          This allows to ignore harmonic cavities.
        divider: number of segments to search for ts
        ts_tol: relative tolerance for ts calculation, default is 1.0-9
          If the search fails to find the synchronous phase the tolerance
          can be increased using the keyword argument or by setting
          `at.DConstant.TStol`
    Returns:
        timelag (float):    (ncav,) array of *Timelag* values
        ts (float):         Time difference with the present value
    """

    def singlev(values):
        vals = np.unique(values)
        if len(vals) > 1:
            raise AtError("values not equal for all cavities")
        return vals[0]

    def eq(x, freq, rfv, tl0, u0):
        omf = 2 * np.pi * freq / clight
        if u0 > 0.0:
            eq1 = (np.sum(-rfv * np.sin(omf * (x - tl0))) - u0) / u0
        else:
            eq1 = np.sum(-rfv * np.sin(omf * (x - tl0)))
        eq2 = np.sum(-omf * rfv * np.cos(omf * (x - tl0)))
        if eq2 > 0:
            return np.sqrt(eq1**2 + eq2**2)
        else:
            return abs(eq1)

    if cavpts is None:
        cavpts = get_bool_index(ring, RFCavity)
    if ts_tol is None:
        ts_tol = DConstant.TStol
    u0 = get_energy_loss(ring, method=method) / ring.periodicity
    freq = np.array([cav.Frequency for cav in ring.select(cavpts)])
    rfv = np.array([cav.Voltage for cav in ring.select(cavpts)])
    tl0 = np.array([cav.TimeLag for cav in ring.select(cavpts)])
    try:
        frf = singlev(freq)
        tml = singlev(tl0)
    except AtError:
        ctmax = clight / np.amin(freq) / 2
        tt0 = tl0[np.argmin(freq)]
        bounds = (-ctmax, ctmax)
        args = (freq, rfv, tl0, u0)
        r = []
        for i in range(divider):
            fact = (i + 1) / divider
            r.append(
                least_squares(
                    eq, bounds[0] * fact + tt0, args=args, bounds=bounds + tt0
                )
            )
            r.append(
                least_squares(
                    eq, bounds[1] * fact + tt0, args=args, bounds=bounds + tt0
                )
            )
        res = np.array([ri.fun[0] for ri in r])
        ok = res < ts_tol
        vals = np.array([abs(ri.x[0]).round(decimals=6) for ri in r])
        if not np.any(ok):
            raise AtError("No solution found for Phis: check RF settings") from None
        if len(np.unique(vals[ok])) > 1:
            warn(
                AtWarning("More than one solution found for Phis: check RF settings"),
                stacklevel=2,
            )
        ts = -r[np.argmin(res)].x[0]
        timelag = ts + tl0
    else:
        vrf = np.sum(rfv)
        if u0 > vrf:
            v1 = ring.periodicity * vrf
            v2 = ring.periodicity * u0
            raise AtError(
                f"The RF voltage ({v1:.3e} eV) is lower than "
                f"the radiation losses ({v2:.3e} eV)."
            )
        timelag = clight / (2 * np.pi * frf) * np.arcsin(u0 / vrf)
        ts = timelag - tml
        timelag *= np.ones(ring.refcount(cavpts))
    return timelag, ts


def set_cavity_phase(
    ring: Lattice,
    *,
    method: ELossMethod = ELossMethod.TRACKING,
    refpts: Refpts | None = None,
    cavpts: Refpts | None = None,
    copy: bool = False,
) -> None:
    """
    Adjust the TimeLag attribute of RF cavities based on frequency,
    voltage and energy loss per turn, so that the synchronous phase is zero.
    An error occurs if all cavities do not have the same frequency.

    .. Warning::

       This function changes the time reference, this should be avoided

    Parameters:
        ring:       Lattice description
        method:     Method for energy loss computation.
          See :py:class:`ELossMethod`.
        cavpts:     Cavity location. If None, use all cavities.
          This allows to ignore harmonic cavities.
        copy:       If True, returns a shallow copy of ring with new
          cavity elements. Otherwise, modify ring in-place.
    """
    # refpts is kept for backward compatibility
    if cavpts is None and refpts is not None:
        warn(FutureWarning('You should use "cavpts" instead of "refpts"'), stacklevel=2)
        cavpts = refpts
    elif cavpts is None:
        cavpts = get_bool_index(ring, RFCavity)
    timelag, _ = get_timelag_fromU0(ring, method=method, cavpts=cavpts)
    set_value_refpts(ring, cavpts, "TimeLag", timelag, copy=copy)


Lattice.get_energy_loss = get_energy_loss
Lattice.energy_loss = property(get_energy_loss)
Lattice.set_cavity_phase = set_cavity_phase
