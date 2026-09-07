# Exact misalignment and curved-frame corrections for the `Exact*` pass methods

Working notes for a possible PR against
[atcollab/at#1067](https://github.com/atcollab/at/pull/1067) (`Quad_curvature_correction`).
Base commit: `51ba4327`. Branch: `xsuite-exact-misalign`.

**Status: not submitted.** Some items below are validated to machine precision,
others are open. The open ones are marked as such and should not be presented
as finished work.

---

## 1. Why

AT applies element misalignments through the linearised `R1`/`R2` matrices and
`T1`/`T2` vectors built by `at.lattice.transformation.transform_elem`. The
*geometry* those encode is exact; their *action on the particle* is truncated to
first order in the transverse momenta. That is fine for the expanded pass
methods, but it is inconsistent with an exact integrator, and it is the bulk of
the AT-vs-Xsuite disagreement for misaligned exact elements.

Measured on a misaligned sector bend, AT's linear treatment disagrees with
Xsuite by **1e-4 to 4e-3** — dominated by `tilt` and `yaw`, and growing with the
rotation angle. That is far above the ~2e-11 integration noise floor.

## 2. What changed

### 2.1 `atintegrators/exact_misalign.h` (new)

Exact rigid-body 6-DOF misalignment, a port of Xsuite's
`track_misalignments.h` into AT's coordinate convention. Wired into
`ExactMultipolePass`, `ExactSectorBendPass` and `ExactRectangularBendPass`.

Design points:

* It reads AT's **existing** attributes — `dx, dy, dz, tilt, pitch, yaw,
  tilt_frame` — rather than introducing new ones. These are plain float
  properties and `PyObject_GetAttrString` resolves properties, so
  `atGetOptionalDouble` reaches them unchanged.
* Opt-in per element via `ExactMisalign=1`, which also **skips `R1`/`T1`**
  (otherwise the misalignment would be applied twice). Default off, so no
  existing result moves.
* The anchor is `Length/2` for `ReferencePoint.CENTRE` (AT's default) and `0`
  for `ENTRANCE`. `ReferencePoint` is a Python `Enum` and is *not* readable
  from C as a number, so an optional `MisalignAnchor` double is exposed for
  the rare non-default case.

Conventions the port depends on (all verified, not assumed):

| | |
|---|---|
| `theta` (Xsuite `rot_y_rad`) | `yaw` |
| `phi` (Xsuite `rot_x_rad`) | **`-pitch`** (the sign flip `at.load.xsuite` already uses) |
| `psi` (Xsuite `rot_s_rad_no_frame`) | `tilt` |
| AT `Yrot(phi)` | Xsuite `YRotation(-phi)` (same reversal for x) |
| AT `r6[ct_]` | absolute path length; Xsuite advances `zeta`, so every longitudinal increment is negated |

AT's arc convention already matches Xsuite's: `transform_elem` places the exit
at `OPp = ((cos a - 1)/h, 0, sin a/h)`, identical to Xsuite's
`matrix_first_part` columns. No arc-sense flip is needed.

### 2.2 `atintegrators/kick_k1h_kn.h` — K0h curvature correction

In the curved frame the multipole potential carries the metric factor
`(1 + irho*x)`, giving one correction per field order. The K1h term was
present; the corresponding dipole term was not, so a `PolynomB[0]` error (or a
`KickAngle` corrector) was integrated without its curvature correction:

```c
r6[1] -= L * irho * (B[0] + B0) * x;   /* H = 1/2 * irho * b0 * x^2, MAD8 eq. 5.15 */
```

`B0` is the `KickAngle`-derived component, so `B[0] + B0` is the total dipole
*error*; the design bend is handled by the exact bend propagator and must not
be included. This applies **only** to `ExactSectorBendPass` — the other two
integrate in a straight frame via `strthinkick`, where no metric factor arises.

### 2.3 `atintegrators/ExactSectorBendPass.c` — dipole fringe sees the total field

Xsuite's `track_magnet_edge.h` builds the field it hands to the dipole fringe as

```c
k0 = knorm[0] + knl[0]/length;    // design + error
```

whereas AT passed `irho` (design only), so `PolynomB[0]` never reached
`bend_fringe`. Changed to pass `irho + B[0] + B0`. Note the *wedge* uses the
design field alone in both codes (`Wedge_single_particle(..., knorm[0])`), so
`bend_edge()` keeps `irho`; and AT's `multipole_fringe` is already called with
`skip_b0=1`, matching Xsuite's `min_order=1`.

### 2.4 Bug fixes found along the way

| file | bug |
|---|---|
| `ExactSectorBendPass.c` | the **exit** quadrupole wedge tested `FringeQuadEntrance` instead of `FringeQuadExit` |
| `pyat/at/lattice/elements/rectangular_bend.py:55` | `float(fsolve(...))` on a 1-element array — NumPy >= 2.0 removed that implicit conversion, so **`rbendtune()` raises `TypeError` on any current NumPy** and a rectangular bend carrying multipoles cannot be set up at all. Fixed with `fsolve(...)[0]`. |

## 3. Results

Harness: `xs_compare.py`. The AT lattice is the single source of truth — it is
converted with `at.line_from_lattice(match_model=True)` and both codes track the
same 11x11 grid at 5 mm, so a mismatch is a real disagreement rather than a
harness artefact. Reference floor for a clean element: **2.0e-11** (sector bend,
`NumIntSteps=200`), **4.9e-14** (multipole).

### 3.1 Misalignment (sector bend, all six DOFs)

| DOF | linear `R1`/`T1` | exact | gain |
|---|---|---|---|
| dx 1 mm | 3.163e-08 | 2.004e-11 | 1.6e3 |
| dy 1 mm | 2.406e-11 | 2.406e-11 | — (see note) |
| dz 1 mm | 5.499e-08 | 2.006e-11 | 2.7e3 |
| tilt 1 mrad | 4.536e-04 | 2.093e-11 | 2.2e7 |
| pitch 1 mrad | 1.254e-07 | 2.054e-11 | 6.1e3 |
| yaw 1 mrad | 1.369e-04 | 2.006e-11 | 6.8e6 |
| tilt 10 mrad | 4.536e-03 | 2.880e-11 | 1.6e8 |
| all six | 4.571e-04 | 2.543e-11 | 1.8e7 |
| all six, large | 4.463e-03 | 2.213e-11 | 2.0e8 |

Note: `dy` on a horizontal bend commutes with the arc, so the exact and linear
transforms are *identically* equal there. It cannot discriminate, and its "1x"
is expected rather than a failure.

### 3.2 Field errors, per order, alone and combined with misalignment

Sector bend, error `1e-4` relative, misalignment = all six DOFs (large set):

| error | alone | + misaligned |
|---|---|---|
| *(misalignment only)* | — | 2.213e-11 |
| B1, A1, B2, A2, B3, A3 | 2.005e-11 | **2.213e-11** |
| B0 | 1.522e-06 | 2.563e-06 |
| A0 | 2.525e-07 | 6.031e-07 |

Orders 1-3 land exactly on the misalignment-only value: **no field x
misalignment cross-terms**. `ExactMultipolePass` is clean at every order
including A0/B0, alone and misaligned (~5.9e-14) — a straight element has no
curved-frame edge.

### 3.3 B0 root cause

The residual is exactly linear in B0 (`1e-5`->1.521e-07, `1e-4`->1.522e-06,
`1e-3`->1.523e-05), sits in `y`/`py`, and **vanishes at face angle 0 with
fringes off** (3.185e-11 = the floor). So the curved-frame body kick, including
the new K0h term, is correct; the disagreement was in the edge. See 2.3.

## 4. Open items — do not present these as done

* **A0 (skew dipole) is not root-caused.** The face/fringe isolation that
  identified the B0 mechanism has not been repeated for A0. The standing
  hypothesis — that Xsuite's `Bend` never computes a skew-dipole fringe at all
  (`MultFringe` is called with `min_order=1`, `DipoleFringe` takes only a scalar
  `k0`, and `bend.h` hardcodes `k0s=0`), making this a physics-completeness
  difference rather than an AT bug — is **unverified in this round**.
* **The rectangular bend is unmeasured.** Every earlier attempt was invalid
  because `rbendtune()` was crashing (2.4). Re-running.
* **The 2.3 fringe fix is not yet validated.** Numbers pending.
* The `MisalignAnchor` question for rectangular bends is unresolved: Xsuite's
  `RBend` anchors on `length_straight` (the chord) while the AT side currently
  defaults to `Length/2` (the arc).
* `ExactPitch`/`ExactYaw` sign conventions are exercised only through
  `transform_elem`; a lattice that sets `R1`/`T1` directly is untested.
* The MEX entry points were updated for the new signature but not exercised.

## 5. Converter bugs (separate from this branch's C changes)

These are in `pyat/at/load/xsuite.py`, are independent of everything above, and
are arguably more important for users than the misalignment work, because they
bite on **default** settings. Confirmed by measurement, on a clean element with
no field errors and no misalignment (face angle 0.5):

| config | AT vs Xsuite |
|---|---|
| `FringeBend=0, FringeQuad=0` | **1.428e-03** |
| `FringeBend=0, FringeQuad=1` | **8.983e-03** |
| face 0, all fringes off | 3.185e-11 (clean) |

1. **The wedge is gated on the wrong flag.** `Bend._set_xs_fringe` derives
   Xsuite's `edge_entry_model` from `FringeQuadEntrance` ("dipole-only" when 0),
   but AT's geometric wedge fires whenever `EntranceAngle != 0`, independent of
   that flag. **AT's default is `FringeQuad*=0`**, so an ordinary dipole with a
   pole-face angle converts to a line that silently disagrees at ~1e-3.
2. **`FringeBendEntrance`/`FringeBendExit` are not mapped at all** — `grep
   FringeBend pyat/at/load/xsuite.py` returns nothing — so switching AT's bend
   fringe off leaves Xsuite's on.

## 6. Caveats on the measurements

* The test dipoles force `FringeQuadEntrance = FringeQuadExit = 1`, which is
  **not** AT's default (0). This was deliberate, to exercise the full edge
  machinery, but it means the tables above do not characterise the
  out-of-the-box dipole — and the default configuration is precisely the one
  that trips converter bug 5.1.
* `FringeQuad*` is not a free knob: the converter derives Xsuite's edge model
  from it, so toggling it moves *both* codes. `compare(..., ref_elem=)` pins the
  Xsuite side and varies only AT when one code must be isolated.
* The 2e-11 floor is set by the drift-kick-drift integration at
  `NumIntSteps=200`, not by the misalignment, so these tests cannot resolve a
  misalignment error below ~1e-11.
* The AT `ct` <-> Xsuite `zeta` mapping used by the harness is a plain sign
  flip, exact only at `delta=0` (where `rvv=1`). All tests run on-momentum.

## 7. Files

| file | role |
|---|---|
| `atintegrators/exact_misalign.h` | the transform (new) |
| `atintegrators/kick_k1h_kn.h` | K0h correction |
| `atintegrators/Exact{Multipole,SectorBend,RectangularBend}Pass.c` | wiring, fringe field, wedge typo |
| `pyat/at/lattice/elements/rectangular_bend.py` | `rbendtune` NumPy 2 fix |
| `pyat/test/test_exact_misalign.py` | regression tests, skipped without xtrack |
| `xs_compare.py`, `t0*.py` | scratch harness — **not for upstream** |

Before submitting, the scratch harness must be separated from the passmethod
changes, and the converter fixes (section 5) should probably be their own PR
since they are independent of the exact-misalignment work.
