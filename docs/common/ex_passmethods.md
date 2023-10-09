# "Exact" passmethods

The "Exact" passmethods are based on E. Forest's book[^Forest]. They are similar to the
tracking in PTC and are not restricted by the small angle approximation.
They use the full expansion of the longitudinal momentum:
```{math}

p_z=\sqrt{(1+\delta)^2-p_x^2-p_y^2}
```
Resulting in:
```{math}

x' &= \frac{p_x}{\sqrt{(1+\delta)^2-p_x^2-p_y^2}} \\
y' &= \frac{p_y}{\sqrt{(1+\delta)^2-p_x^2-p_y^2}}
```
"Exact" passmethods are more computationally intensive and much slower than
the default methods. They are mainly useful for small rings (small bending
radius, large angles).

## `ExactDriftPass`
Exact integration in a free space region

(ExactMultipolePass)=
## `ExactMultipolePass`

(ExactSectorBendPass)=
## `ExactSectorBendPass`
Exact integration in a curved magnet.

This method uses the bend-kick split of the Hamiltonian. The "bend" step integrates
the order 0 of the field expansion (dipole field) while the kick includes the effects
of the higher orders of the field expansion and of the synchrotron radiation.
Following the notations in [^Forest],
the map corresponds to {math}`\mathcal{Y}(\varepsilon_1)
\mathcal{F}_1\mathcal{U}(-\varepsilon_1)\mathcal{W}\mathcal{U}(-\varepsilon_2)
\mathcal{F}_2\mathcal{Y}(\varepsilon_2)` with:
- {math}`\mathcal{Y}(\varepsilon_1)`: y-axis rotation (Eq. 10.26)
- {math}`\mathcal{F}_1`: dipole fringe field in the hard-edge limit (Eq. 13.13), [^F2]
- {math}`\mathcal{U}(-\varepsilon_1)`: entrance wedge (Eq. 12.41)
- {math}`\mathcal{W}`: bend-kick sequence in cylindrical geometry(Eq. 12.18)
- {math}`\mathcal{U}(-\varepsilon_2)`: exit wedge
- {math}`\mathcal{F}_2`: dipole fringe field in the hard-edge limit
- {math}`\mathcal{Y}(\varepsilon_2)`: y-axis rotation

```{Tip}
For a pure dipole field, without synchrotron radiation, it is possible to
integrate the whole magnet in one step, resulting in a much faster transfer. For
this, use `NumIntSteps=0`.
```

PolynomB, PolynomA, MaxOrder, NumIntSteps
: see [StrMPoleSymplectic4Pass](#StrMPoleSymplectic4Pass)

Length, BendingAngle, EntranceAngle, ExitAngle
: see [BndMPoleSymplectic4Pass](#BndMPoleSymplectic4Pass)

```{Caution}
`ExactSectorBendPass` show a small discontinuity around origin. Therefore it is
not recommended to use it for computations based on transfer matrices (linear
optics). This can be mitigated by increasing the differentiation steps `XYStep`
and `DPStep` with respect to their default values.
```

(ExactRectangularBendPass)=
## `ExactRectangularBendPass`
Exact integration in a bending magnet with Cartesian layout.

This method uses the drift-kick split of the Hamiltonian in the Cartesian
coordinates of the magnet.
Following the notations in [^Forest],
the map corresponds to {math}`\mathcal{Y}(\varepsilon_1)
\mathcal{F}_1\mathcal{U}(\theta/2-\varepsilon_1)\mathcal{D}
\mathcal{U}(\theta/2-\varepsilon_2)\mathcal{F}_2\mathcal{Y}(\varepsilon_2)` with:
- {math}`\mathcal{Y}(\varepsilon_1)`: y-axis rotation (Eq. 10.26)
- {math}`\mathcal{F}_1`: dipole fringe field in the hard-edge limit (Eq. 13.13), [^F2]
- {math}`\mathcal{U}(\theta/2-\varepsilon_1)`: entrance wedge (Eq. 12.41)
- {math}`\mathcal{D}`: drift-kick sequence
- {math}`\mathcal{U}(\theta/2-\varepsilon_2)`: exit wedge
- {math}`\mathcal{F}_2`: dipole fringe field in the hard-edge limit
- {math}`\mathcal{Y}(\varepsilon_2)`: y-axis rotation

For consistency with other passmethods, the `Length` attribute is the length {math}`L`
of the arc within the magnet, and not the Cartesian length {math}`L_c` of the
magnet. {math}`L=L_c\frac{\theta/2}{sin(\theta/2)}`

If the magnet field includes quadrupole or higher components, the reference trajectory 
in the magnet is no more an arc of a circle. A tuning of the `X0ref` attribute is
necessary to get the correct output angle from the magnet. This can be seen as a
horizontal translation of the magnet until the correct deviation angle is obtained.
Similarly, the path lengthening must be adjusted through the `RefDZ` attribute to
take into account the length of the non-circular trajectory.

This tuning is performed using a dedicated function/method:

- python: {py:meth}`.Dipole.rbendtune`
  ```python
  # Identify the rectangular bends (for example...)
  rbends = ring.get_bool_index(checkattr("PassMethod", "ExactRectangularBendPass")
  # Set their correct attributes
  for dip in ring.select(rbends):
      dip.rbendtune()
  ```
- Matlab: `atrbendtune()`
  ```Matlab
  % Identify the rectangular bends (for example...)
  rbends=atgetcells(ring,'PassMethod', 'ExactRectangularBendPass');
  % Set their correct attributes
  ring(rbends)=cellfun(@attunerbend,ring(rbends),'UniformOutput',false);
  ```

PolynomB, PolynomA, MaxOrder, NumIntSteps
: see [StrMPoleSymplectic4Pass](#StrMPoleSymplectic4Pass)

Length, BendingAngle, EntranceAngle, ExitAngle
: see [BndMPoleSymplectic4Pass](#BndMPoleSymplectic4Pass)

[^Forest]: Étienne Forest, _Beam Dynamics, a new Attitude and Framework_, 
Harwood Academic Publishers.

[^F2]: É. Forest, S.C. Leemann, F Schmidt, _Fringe Effects in MAD - Part I_.