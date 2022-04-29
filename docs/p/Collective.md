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


### Wake element

The `WakeElement` is the **AT** element to be integrated in the lattice, it is using the passmethod `WakefieldPass.c` as the integrator and is initialized with the following arguments:

`args`:
- family name: string name of the element
- ring: lattice object
- wake: wake object

`kwargs`: 
- PassMethod = 'WakeFieldPass'
- NumParticles = 0,   number of particles
- Nslice = 101,  number of slices per bunch
- Nturns = 1,    number of turn for the multi-turn wake field
- ZCuts = None,  limits for fixed slicing, default is adaptive
- NormFact = [1,1,1], normalization for the 3 planes, to account for beta function at the observation point for example

The public attributes/properties for this element are:
- NumParticles, number of charges
- Normfact, normalization factor
- Zcuts, fixed slicing limits
- WakeT/DX/DY/QY/QY/Z: component of the wake considered for tracking, None is returned in case they are defined
- Nturns, number of turns considered for the multi-turn wake, starts at 1 (multi-turn disabled)
- Current, beam current [A]

Functions are available to re-initialize some attributes:
- `WakeElement.rebuild_wake(wake)`: rebuilds the wake field using a new `Wake` object
- `WakeElement.clear_history()`: clear the turn history used in the multi-turn wake calculation
- `WakeElement.set_normfactxy()`: applies a normalization by the beta-function at the observation point

Similarly to the `Wake` object specific cases are provided: `ResonatorElement`, `LongResonatorElement`, `ResWallElement`. In these case the `Wake` Object is constructed internally and the proper `args` and `kwargs` to build it need to be provided in the initialization of the `WakeElement`. These arguments are defined as properties and can be accessed and changed by the user after initialization.

### Example usage

Practical examples are found in `pyat/examples/CollectiveEffects`. The following shows an example of tracking with a transverse resonator. It uses the lattice `esrf.m` that can be found in `at/machine_data`.

The first step is to load the lattice, set the rf frequency to its nominal value, turn radiations off and generate a `fast_ring`. The `fast_ring` reduces the lattice to a few objects: a 6x6 transfer map, an RF cavity, an element to model chromaticity and amplitude detuning and a diffusion element in case the radiations are turned on.

```python
import at
from at.collective import Wake, ResonatorElement, WakeComponent

ring = at.load_m('esrf.m')
ring.radiation_off(cavity_pass='RFCavityPass')
ring.set_rf_frequency()
fring, _ = at.fast_ring(ring)
```

The `WakeElement` is then created and added to the `fast_ring`. In this simulation we consider only single turn effects, `srange` is therefore build accordingly, for single turn short range wakes the last 3 arguments are ignored, the interpolation table is defined from 0 to 0.3m in steps of 1.0e-5m. A 1 GHz vertical broadband resonator is then created, a current of 10mA is assigned to the element and it is appended to the lattice.

```python
srange = Wake.build_srange(0., 0.3, 1.0e-5, 1, ring.circumference, ring.circumference)

fr = 1.0e9
qfactor = 1
rshunt = 6e6
welem = ResonatorElement('wake', ring, srange, WakeComponent.DY, fr, qfactor, rshunt)
welem.Current = 10e-3

fring.append(welem)
```

A beam of 10000 macro-particles is then generated and track through the lattice for 1000 turns

```python
sigm = at.sigma_matrix(ring.radiation_on(copy=True))
part = at.beam(10000, sigm)
part_out = at.lattice_pass(fring, part, nturns=1000)
```
**Warning:** AT will save the 6 coordinates of all particles at each turn in part_out, the `C` tracking engine therefore needs to allocate the associated memory, for very large arrays this may cause allocation errors because of insufficient memory. In order to prevent such error it is possible to track one turn at a time and save intermediate values at each turn

```python
xmean = numpy.zeros(1000)
for i in range(1000):
    _ = at.lattice_pass(fring, part, nturns=1)
    x_mean[i] = numpy.mean(part[0, :], axis=0)
```
