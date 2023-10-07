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
These methods are fast and the default methods assigned when creating AT elements.

## `DriftPass`
Field-free space

```{rubric} Default passmethods
```

## `StrMPoleSymplectic4Pass`
4th order Forest-Knuth integrator for straight magnets

Length
: Length of the arc.

PolynomB, PolynomA
: Polynomial field expansion.

MaxOrder
: Maximum order of the polynomial expansion: only indices from 0 to MaxOrder
(included) are used.

NumIntSteps
: Number of integration steps (magnet slices). Optional,default 10.

## `BndMPoleSymplectic4Pass`
4th order Forest-Knuth integrator for curved magnets

Length
: Length of the arc

BendingAngle
: Dipole bending angle

PolynomB, PolynomA
: Polynomial field expansion 

```{rubric} Linear passmethods
```

## `QuadLinearPass`

## `BendLinearPass`
