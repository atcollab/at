"""Tracking-fidelity harness: one AT element vs its Xsuite conversion.

The AT lattice is the single source of truth: it is converted with
``at.line_from_lattice`` and both codes track the same particle grid.  No
misalignment is applied by hand on the Xsuite side -- whatever the converter
produces is what gets compared, so a mismatch is a real AT/Xsuite disagreement
(or a converter bug) rather than an artefact of the harness.
"""

from __future__ import annotations

import numpy as np
import at
import xobjects
import xtrack as xt

# No prebuilt kernels in this venv: allow just-in-time compilation.
xobjects.context_cpu.allow_no_prebuilt_kernel = True

ENERGY = 10e9
NAMES = ["x", "px", "y", "py", "delta", "ct"]


def make_particles(n: int, amp: float, planes=(0, 2), energy=ENERGY, delta=0.0):
    """Grid of n x n particles in the two given planes."""
    a = np.linspace(-amp, amp, n)
    g1, g2 = np.meshgrid(a, a)
    parts = np.zeros((6, g1.size))
    parts[planes[0]] = g1.ravel()
    parts[planes[1]] = g2.ravel()
    parts[4] = delta
    return parts


def track_at(lattice, parts):
    pin = np.asfortranarray(parts.copy())
    pout, *_ = lattice.track(pin)
    # (6, n_particles, n_refpts, n_turns) -> (6, n_particles)
    return np.asarray(pout).reshape(6, parts.shape[1], -1)[:, :, -1]


def track_xs(line, parts, energy=ENERGY):
    if line.particle_ref is None:
        line.particle_ref = xt.Particles(p0c=energy, mass0=xt.ELECTRON_MASS_EV)
    p = line.build_particles(
        x=parts[0], px=parts[1], y=parts[2], py=parts[3],
        delta=parts[4], zeta=parts[5],
    )
    line.track(p)
    out = np.full((6, parts.shape[1]), np.nan)
    idx = p.particle_id[p.state > 0]
    out[0, idx] = p.x[p.state > 0]
    out[1, idx] = p.px[p.state > 0]
    out[2, idx] = p.y[p.state > 0]
    out[3, idx] = p.py[p.state > 0]
    out[4, idx] = p.delta[p.state > 0]
    # AT stores path lengthening, Xsuite the (oppositely signed) zeta
    out[5, idx] = -p.zeta[p.state > 0]
    return out


def build_line(lattice):
    line = at.line_from_lattice(lattice, match_model=True)
    if line.particle_ref is None:
        line.particle_ref = xt.Particles(p0c=lattice.energy,
                                         mass0=xt.ELECTRON_MASS_EV)
    return line


def compare(elem, amp=5e-3, n=11, planes=(0, 2), energy=ENERGY, delta=0.0,
            label="", verbose=True, show_xs=False):
    """Track one element through both codes; return the per-coordinate max |diff|."""
    lattice = at.Lattice([elem.copy()], energy=energy, periodicity=1)
    line = build_line(lattice)

    if show_xs:
        e0 = line[0]
        fields = ["shift_x", "shift_y", "shift_s", "rot_s_rad_no_frame",
                  "rot_x_rad", "rot_y_rad", "rot_shift_anchor",
                  "edge_entry_model", "edge_exit_model"]
        got = {f: getattr(e0, f, None) for f in fields}
        print("   xsuite:", {k: v for k, v in got.items() if v})

    parts = make_particles(n, amp, planes, energy, delta)
    a_out = track_at(lattice, parts)
    x_out = track_xs(line, parts, energy)

    diff = a_out - x_out
    finite = np.isfinite(diff).all(axis=0)
    if not finite.any():
        print(f"   {label}: ALL PARTICLES LOST")
        return None
    dmax = np.nanmax(np.abs(diff[:, finite]), axis=1)
    if verbose:
        worst = np.max(dmax)
        print(f"   {label}: max|AT-Xsuite| = {worst:.3e}   "
              f"({finite.sum()}/{len(finite)} kept)")
        print("      " + "  ".join(f"{n}={v:.2e}" for n, v in zip(NAMES, dmax)))
    return dmax


def sector_bend(length=2.0, angle=1.0, k1=1.0, nsteps=200, face=0.5,
                fringe_bend=1, fringe_quad=1, **kwargs):
    d = at.Dipole("D", length, angle, k=k1, PassMethod="ExactSectorBendPass",
                  NumIntSteps=nsteps, EntranceAngle=face, ExitAngle=face,
                  FringeBendEntrance=fringe_bend, FringeBendExit=fringe_bend,
                  FringeQuadEntrance=fringe_quad, FringeQuadExit=fringe_quad)
    for k, v in kwargs.items():
        setattr(d, k, v)
    return d


def multipole(length=1.0, k1=0.5, nsteps=200, **kwargs):
    d = at.Multipole("M", length, [0.0, 0.0, 0.0, 0.0], [0.0, k1, 0.0, 0.0],
                     PassMethod="ExactMultipolePass", NumIntSteps=nsteps)
    for k, v in kwargs.items():
        setattr(d, k, v)
    return d


def rect_bend(length=2.0, angle=0.4, k1=0.3, nsteps=200, **kwargs):
    d = at.Dipole("R", length, angle, k=k1,
                  PassMethod="ExactRectangularBendPass", NumIntSteps=nsteps,
                  EntranceAngle=angle / 2, ExitAngle=angle / 2,
                  FringeBendEntrance=1, FringeBendExit=1,
                  FringeQuadEntrance=1, FringeQuadExit=1)
    for k, v in kwargs.items():
        setattr(d, k, v)
    return d


def apply_misalign(elem, exact=True, **kwargs):
    """Apply an AT misalignment and (optionally) request the exact treatment."""
    at.transform_elem(elem, **kwargs)
    if exact:
        elem.ExactMisalign = 1
    return elem
