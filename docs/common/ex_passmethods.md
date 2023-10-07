# "Exact" passmethods

The "Exact" passmethods are base on E. Forest's book _"Beam Dynamics,
a new Attitude and Framework"_. They are similar to the tracking in PTC
and are not restricted by the small angle approximation.
They use the full expansion of the longitudinal momentum:
```{math}

p_z=\sqrt{(1+\delta)^2-p_x^2-p_y^2}
```
Resulting in:
```{math}

x' &= \frac{p_x}{\sqrt{(1+\delta)^2-p_x^2-p_y^2}} \\
y' &= \frac{p_y}{\sqrt{(1+\delta)^2-p_x^2-p_y^2}}
```

## `ExactDriftPass`

## `ExactSectorBendPass`
Exact integration in a curved magnet.

This method uses the drift-kick split of the Hamiltonian. The "bend" step integrates
the order 0 of the field expansion (dipole field) while the kick includes the effect
of the higher orders of the field expansion and of the synchrotron radiation.

For a pure dipole field, without synchrotron radiation, it is therefore possible to
integrate the whole magnet in one step.

## `ExactRectangularBendPass`
Exact integration in a bending magnet with Cartesian layout.

For consistency with other passmethods, the Length attribute is the length of the arc
within the magnet, and not the Cartesian length of the magnet.

If the magnet field includes quadrupole or higher components, the reference trajectory 
in the magnet is no more an arc of a circle. A tuning of the `X0ref` attribute is
necessary to get the correct output angle from the magnet. This can be seen as a
horizontal translation of the magnet until the correct deviation angle is obtained.

This tuning is performed using a dedicated function/method:

- python
  ```{code} python
  # Identify the rectangular bends (for example...)
  rbends = ring.get_bool_index(checkattr("PassMethod", "ExactRectangularBendPass")
  # Set their correct attributes
  for dip in ring.select(rbends):
      dip.rbendtune()
  ```
- Matlab
  ```{code} Matlab
  % Identify the rectangular bends (for example...)
  rbends=atgetcells(ring,'PassMethod', 'ExactRectangularBendPass');
  % Set their correct attributes
  ring(rbends)=cellfun(@attunerbend,ring(rbends),'UniformOutput',false);
  ```