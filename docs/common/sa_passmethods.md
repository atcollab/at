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

(GWigSymplecticPass)=
## `GWigSymplecticPass`
2 {sup}`nd` and 4{sup}`th` order Forest-Knuth integrator for wigglers without radiation

Length
: Length {math}`L` of the element.

MaxOrder, NumIntSteps
: see [StrMPoleSymplectic4Pass](#StrMPoleSymplectic4Pass)

Period
: Wiggler period {math}`L_w`.

Peak magnetic field, B_0
: Maximum magnetic field {math}`B_0` of the wiggler.

Integration method, Nmeth
: Nmeth indicates the integration method: 2nd or 4th order integrator

NHharm
: Number of horizontal harmonics of the wiggler

NVharm
: Number of vertical harmonics of the wiggler

By
: 6 x NHharm array containing the following quantities. row 1: Number of current horizontal
harmonic; row 2: relative amplitudes of wiggler harmonics; row 3: {math}`k_x*L_w/2\pi`;
row 4: {math}`k_y*L_w/2\pi`; row 5: {math}`k_z*L_w/2\pi`; row 6: {math}`\theta_n`, the relative 
phase of the n{sup}`th` wiggler harmonic. {math}`k_x`, {math}`k_y` and {math}`k_z` are defined in [^Wu].

Bx
: 6 x NVharm array containing the following quantities. row 1: Number of current vertical
harmonic; row 2: relative amplitudes of wiggler harmonics; row 3: {math}`k_x*L_w/2\pi`;
row 4: {math}`k_y*L_w/2\pi`; row 5: {math}`k_z*L_w/2\pi`; row 6: {math}`\theta_n`, the relative 
phase of the n{sup}`th` wiggler harmonic. {math}`k_x`, {math}`k_y` and {math}`k_z` are defined in [^Wu].


[^Wu]: Y. K. Wu, E. Forest, and D. S. Robin, _Explicit symplectic integrator for s-dependent static magnetic
field_, Phys. Rev. E, 68:046502, Oct 2003.