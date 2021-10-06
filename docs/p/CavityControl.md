---
title: Controlling the RF cavities
permalink: /pcavitycontrol.html
---
A lattice may contain multiple RF cavities, grouped according to different RF systems:
main cavities, harmonic cavities…

AT provides simple tools to tune them, with methods and properties of the
`Lattice` object.

### Lattice methods for cavity control

The values handled by these methods concern the full ring (`periodicity` $$\times$$ superperiod).

All methods have a `cavpts` argument, used to select the cavities concerned by the command.
The `cavpts` argument is used as follows:
- `cavpts` is a "refpts" type (integer, integer or boolean array, callable): it is used to select the cavities,
- `cavpts` is `None` (default value), and the `Lattice` object has a `cavpts` attribute: the lattice attribute is used to select the cavities,
- `cavpts` is `None`, and the lattice has no `cavpts` attribute (or it is `None`): all cavities are selected.

So the easier way to use it is:
- **single RF system** (main cavities): forget the `cavpts` argument. The default is to use all cavities,
- **main RF system + harmonic system**: set the `Lattice` `cavpts` attribute to the location of the main cavities,
  so that the default behaviour is to drive the main system. Use the method's `cavpts` argument to drive the harmonic cavities.

All methods also have a `copy` argument to select either in-place modification
of the lattice, or creation of a shallow copy with modified cavities.

#### Voltage:
```voltage = Lattice.get_rf_voltage(cavpts=None)```

```Lattice.set_rf_voltage(voltage, cavpts=None, copy=False)```

The specified voltage is equally distributed among all selected cavities in all superperiods.

#### Frequency:
`frequency = Lattice.get_rf_frequency(cavpts=None)`

`Lattice.get_rf_frequency(frequency=None, dp=None, dct=None, cavpts=None, copy=False)`

If the frequency in `None`, the method will set the frequency to the nominal value,
according to the revolution frequency and harmonic number. An optional
off-momentum may be applied using the `dp` or `dct` arguments. The frequency
shift is then computed using the linear slip factor $$\eta_c = 1/\gamma^2 - \alpha_c$$ ,
so that the resulting `dp` may slighly differ from the specified value.
#### Harmonic number
`harmonic_number = Lattice.get_rf_harmonic_number(cavpts=None)`

`Lattice.set_rf_harmonic_number(harmonic_number, cavpts=None, copy=False)`
#### Time lag
`time_lag = Lattice.get_rf_timelag(cavpts=None)`

`Lattice.set_rf_timelag(time_lag, cavpts=None, copy=False)`

The time lag is expressed in values of path lengthening "c&tau;", the 6<sup>th</sup> particle coordinate [m].
#### All-in-one method
`Lattice.set_cavity(ring, Voltage=None, Frequency=None, HarmNumber=None,
               TimeLag=None, cavpts=None, copy=False)`

This method sets only the explicity provided values, default ones are left unchanged.
For the frequency, a special value `at.Frf.NOMINAL` means nominal frequency,
according to the revolution frequency and harmonic number.

### Lattice properties
The properties provide an even easier way to control the cavities, but are restricted
to the default behaviour of the equivalent Lattice method:
- cavities are selected by the `Lattice` `cavpts` argument (all cavities by default),
- Setting a property modifies the ring in-place (no copy).

`Lattice.rf_voltage`

`Lattice.rf_frequency`

The special value `at.Frf.NOMINAL` means nominal frequency,
according to the revolution frequency and harmonic number.

`Lattice.harmonic_number`

`Lattice.time_lag`