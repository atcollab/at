Cavity Control
==============

A lattice may contain multiple RF cavities, grouped according to 
different RF systems: main cavities, harmonic cavities…
AT provides simple tools to tune them, with methods and properties of
the ``Lattice`` object.

Lattice methods for cavity control
----------------------------------

All methods have ``array`` and ``cavpts`` keyword arguments, used to
select the cavities concerned by the command.

-  if ``array is True``, the output attribute value is an array as long
   as the number of selected cavities. The input argument must be an
   array as long as the number of selected cavities or a scalar which
   will be broadcasted to the number of cavities,
-  if ``array is False`` (default), the input and output are scalars.
   The scalar value applies to the set of cavities with the lowest
   frequency within the selection. The other cavities are ignored in
   ``get_*`` methods. For ``set_*`` methods, the other cavities are
   scaled as explained if the specific method description.

The ``cavpts`` argument is used as follows:

- ``cavpts`` is a “refpts” type (integer, integer or boolean array, callable):
  it is used to select the cavities.
- ``cavpts is None`` (default value), and the ``Lattice`` object has a
  ``cavpts`` attribute: the lattice attribute is used to select the cavities.
- ``cavpts is None``, and the lattice has no ``cavpts`` attribute (or it is
  ``None``): all cavities are taken into account.

.. note::

   -  **single RF system:** you can forget the ``cavpts`` argument: by default
      the methods address all cavities. However, using the ``Lattice.cavpts``
      attribute makes calls significantly faster by skipping the search for
      cavities.
   -  **complex system:** an easy way is to have the ``Lattice.cavpts``
      address the accelerating cavities so that they will be driven by default.
      Harmonic cavities may be driven using the ``cavpts`` argument.

.. note::

   Adding to the Lattice "refpts" type attributes addressing the different
   cavity sets make them available everywhere the lattice is visible.

All ``set_*`` methods also have a ``copy`` argument to select either
in-place modification of the lattice, or creation of a shallow copy with
modified cavities.

Voltage:
~~~~~~~~

``voltage = Lattice.get_rf_voltage(cavpts=None, array=False)``

The scalar voltage is the sum of the cavity voltages of the cavities
with the lowest frequency within the selection, multiplied by the
periodicity.

``Lattice.set_rf_voltage(voltage, cavpts=None, array=False, copy=False)``

For array == False, all the existing voltages are scaled to reach the
specified value on the fundamental mode.

Frequency:
~~~~~~~~~~

``frequency = Lattice.get_rf_frequency(cavpts=None, array=False)``

The frequency of the fundamental mode is returned.

``Lattice.get_rf_frequency(frequency=None, dp=None, dct=None, cavpts=None, array=False, copy=False)``

If the frequency is None, the method will set the frequency to the
nominal value, according to the revolution frequency and harmonic
number. An optional off-momentum may be applied using the ``dp`` or
``dct`` arguments. The frequency shift is then computed using the linear
slip factor :math:`\eta_c = 1/\gamma^2 - \alpha_c`, so that the resulting
``dp`` may slightly differ from the specified value.

For array == False, the value is applied to the fundamental mode
cavities and the frequency of all other cavities is scaled by the same
ratio.

Time lag
~~~~~~~~

The time lag is expressed in values of path lengthening “cτ”, the 6th
particle coordinate [m].

``time_lag = Lattice.get_rf_timelag(cavpts=None, array=False)``

The time lag of the fundamental mode is returned.

``Lattice.set_rf_timelag(time_lag, cavpts=None, array=False, copy=False)``

For array == False, the time lag is applied to the fundamental mode
cavities and the time lag of all the other selected cavities is shifted
by the same amount.

All-in-one method
~~~~~~~~~~~~~~~~~

``Lattice.set_cavity(ring, Voltage=None, Frequency=None, TimeLag=None, cavpts=None, array=False, copy=False)``

This method sets only the explicitly provided values, the other ones are
left unchanged. For the frequency, a special value ``at.Frf.NOMINAL``
means nominal frequency, according to the revolution frequency and
harmonic number.

The behaviour of the ``cavpts`` and ``array`` keywords is the same as
for individual methods.

Lattice properties
------------------

The properties provide an even easier way to control the cavities, but
are restricted to the default behaviour of the equivalent Lattice
method: - cavities are selected by the ``Lattice.cavpts`` attribute (all
cavities by default), - Setting a property modifies the ring in-place
(no copy).

``Lattice.rf_voltage``

``Lattice.rf_frequency``

The special value ``at.Frf.NOMINAL`` means nominal frequency.

``Lattice.harmonic_number``

``Lattice.time_lag``
