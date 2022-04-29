---
title: Collective Effects
---

A collective effects subpackage `pyat/at/collective` allows to model impedance driven collective effects and perform multi-particle tracking. Presently only single bunch effects are included, however multi-turn wake are included to model beam loading effects for example. It is possible to build pyAT with MPI to perform computing intensive simulations on a cluster.
The bunch is sliced longitudinal and the wake field from leading slices is applied on trailing slices iteratively, at each turn the slicing is updated and the kick recomputed.

### Wake object

The `Wake` provides an interface to create an object containing the wake field information. It is then passed to the lattice element that is used for tracking. A `Wake` is defined by its `s` coordinate its wake componenents: transverse dipole, transverse quadrupole and longitudinal and its type: custom, resistive wall or resonator. The type and components are determined by enumerators:

`Wake.WakeComponent`:
- `DX`: horizontal dipole
- `DY`: vertical dipole
- `QX`: horizontal quadrupole (or detuning)
- `QY`: vertical quadrupole (or detuning)
- `Z`: longitudinal

`Wake.WakeType`:
- `FILE`: import the wake from a file
- `TABLE`: provide vectors
- `RESONATOR`: analytical resonator
- `RESWALL`: analytical resistive wall

The wake object is then built using:

```python
wake = Wake(srange)
wake.add(WakeType,WakeComponent, *args, *kwargs)
```
`srange` is user defined and defined the points of the interpolation table. It is nevertheless it is more efficient to provide points only where needed to avoid storing large double array. A static method is provided to build the `srange`:
```python
srange = Wake.build_srange(start, bunch_ext, short_step, long_step, bunch_interval, totallength)
```
where:
- `start`           starting s-coordinate of the table (can be negative for wake potential)
- `bunch_ext`       maximum bunch extension, function generates data at +/- `bunch_ext` around the bucket center
- `short_step`      step size for the short range wake table
- `long_step`       step size for the long range wake table
- `bunch_interval`  minimum bunch interval data will be generate for each bunch_interval step
- `totallength`     total length of the wake table, has to contain the full bunch extension

It is possible to add as many components as needed, all the componenents are resampled (interpolated) on the provided `srange`. `args` and `kwargs` allow to provide inputs specific to each `WakeType` that are defined as follows:

`WakeType.FILE`:
- `args`: file name
- `kwargs`: `scol=0` column for the `s` coordinate, `wcol=1` column for the wake field, `sfact=1` and `wfact=1` column for the wake field scaling factors, `delimiter=None` column delimiter, 'skiprows=0' allows to skip headers. See `numpy.loadtxt()` for details

`WakeType.TABLE`:
- `args`: vector of `s` coordinate, vector of wake field

`WakeType.RESONATOR`:
- `args`: frequency, Q factor, shunt impedance, relativistic beta
- `kwargs`: `yokoya_factor=1'

`WakeType.RESWALL`:
- `args`: length, vacuum pipe radius, conductivity, relativistic beta
- `kwargs`: `yokoya_factor=1'

Static factory methods `Wake.resonator`, `Wake.long_resonator` and `Wake.resistive_wall` are provided for quick instanciating of analytical cases.
These can be called as:
```python
wake = Wake.resonator(WakeComponent, *args, *kwargs)
```
Vector inputs are allows in which case `nelems=1` has to be specified, each inputs is then either a scalar or a vector of shape `(nelems,)`.
