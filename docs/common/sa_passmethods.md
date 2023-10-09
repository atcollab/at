# "Small angle" Passmethods

The "Small angle" passmethods use a linearised approximation of the longitudinal
momentum:
```{math}

p_z=1+\delta
```
resulting in:
```{math}

x' &= \frac{p_x}{1+\delta} \\ y' &= \frac{p_y}{1+\delta}
```
These methods are fast and are assigned by default when creating AT elements.

```{rubric} Default passmethods
```

## `IdentityPass`
Thin element with no effect. However, the element may define transverse apertures
limiting the physical aperture,

## `DriftPass`
Field-free space in the small angle approximation.

Length
: Drift length.

## `RFCavityPass`
RF Cavity

Length
: Cavity length. If `Length` is non-zero, a drift of half-length is added on each
side of a thin cavity.

Voltage
: Cavity voltage

Frequency
: Cavity frequency

TimeLag
: Time lag expressed as path lengthening: {math}`\beta c \tau`.

PhaseLag
: Phase lag in radians. `TimeLag` and `PhaseLag` are accumulated as:
  
  ```{math}
  \mathbf{PhaseLag} - 2 \pi f_{RF} \frac{\mathbf{TimeLag}}{\beta c}
  ```

(StrMPoleSymplectic4Pass)=
## `StrMPoleSymplectic4Pass`
4{sup}`th` order Forest-Knuth integrator for straight magnets

Length
: Magnet Length.

PolynomB, PolynomA
: Polynomial field expansion.

MaxOrder
: Maximum order of the polynomial expansion: only indices from 0 to `MaxOrder`
(included) are used.

NumIntSteps
: Number of integration steps (magnet slices). Optional, default 10.

(BndMPoleSymplectic4Pass)=
## `BndMPoleSymplectic4Pass`
4{sup}`th` order Forest-Knuth integrator for curved magnets

Length
: Length {math}`L` of the arc.

PolynomB, PolynomA, MaxOrder, NumIntSteps
: see [StrMPoleSymplectic4Pass](#StrMPoleSymplectic4Pass)

BendingAngle
: Dipole bending angle {math}`\theta`

EntranceAngle
: Angle of the entrance pole face {math}`\varepsilon_1` with respect to the plane
perpendicular to the input trajectory. Use 0 for a sector magnet,
{math}`\theta/2` for a rectangular magnet.

EntranceAngle
: Angle of the exit pole face {math}`\varepsilon_2` with respect to the plane
perpendicular to the output trajectory. Use 0 for a sector magnet,
{math}`\theta/2` for a rectangular magnet.

## `ThinMPolePass`
Thin multipolar kick

PolynomB, PolynomA, MaxOrder
: see [StrMPoleSymplectic4Pass](#StrMPoleSymplectic4Pass)

```{rubric} Linear passmethods
```
These passmethods are kept for backward compatibility but their use is discouraged.

## `QuadLinearPass`

Length
: Quadrupole Length.

PolynomB, K
: Quadrupole strength. If `PolynomB[1]` exists, it is used for the strength, otherwise
`K` is used

(BendLinearPass)=
## `BendLinearPass`

Length
: Length {math}`L` of the arc.

PolynomB, K
: Quadrupole strength. If `PolynomB[1]` exists, it is used for the strength, otherwise
`K` is used

BendingAngle
: Dipole bending angle {math}`\theta`

EntranceAngle
: Angle of the entrance pole face {math}`\varepsilon_1` with respect to the plane
perpendicular to the input trajectory. Use 0 for a sector magnet,
{math}`\theta/2` for a rectangular magnet.

ExitAngle
: Angle of the exit pole face {math}`\varepsilon_2` with respect to the plane
perpendicular to the output trajectory. Use 0 for a sector magnet,
{math}`\theta/2` for a rectangular magnet.
