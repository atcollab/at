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
{#anchor1}
These methods are fast and the default methods assigned when creating AT elements.

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

(strmpole)=
## `StrMPoleSymplectic4Pass`
4{sup}`th` order Forest-Knuth integrator for straight magnets

Length
: Length {math}`L` of the arc.

PolynomB, PolynomA
: Polynomial field expansion.

MaxOrder
: Maximum order of the polynomial expansion: only indices from 0 to MaxOrder
(included) are used.

NumIntSteps
: Number of integration steps (magnet slices). Optional,default 10.

(bndmpole)=
## `BndMPoleSymplectic4Pass`
4{sup}`th` order Forest-Knuth integrator for curved magnets

Length, PolynomB, PolynomA, MaxOrder, NumIntSteps
: see [here](#strmpole)

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

```{rubric} Linear passmethods
```

## `QuadLinearPass`

## `BendLinearPass`
