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
   will be broadcasted to all the selected cavities,
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

   In most cases, the default behaviour (``cavpts=None``, ``array=False``) gives
   the expected result: the scalar input arguments drives the main cavivities, and
   all the other cavities are scaled or shifted.

   For a more detailed control, the easiest way is to use ``cavpts`` to select
   groups of similar cavities which can be tuned together.

All ``set_*`` methods also have a ``copy`` argument to select either
in-place modification of the lattice, or creation of a shallow copy with
modified cavities.

Voltage
~~~~~~~

:py:func:`voltage = ring.get_rf_voltage(cavpts=None, array=False) <.get_rf_voltage>`

The scalar voltage is the sum of the cavity voltages of the cavities
with the lowest frequency within the selection, multiplied by the
periodicity.

:py:func:`ring.set_rf_voltage(voltage, cavpts=None, array=False, copy=False) <.set_rf_voltage>`

When ``array is False``, the voltages of all selected cavities are scaled to reach
the specified value on the fundamental mode.

Frequency
~~~~~~~~~

:py:func:`frequency = ring.get_rf_frequency(cavpts=None, array=False) <.get_rf_frequency>`

The frequency of the fundamental mode is returned.

:py:func:`ring.set_rf_frequency(frequency=None, dp=None, dct=None, cavpts=None, array=False, copy=False) <.set_rf_frequency>`

If the frequency is None, the method will set the frequency to the
nominal value, according to the revolution frequency and harmonic
number. An optional off-momentum may be applied using the ``dp``, ``dct``
or ``df`` arguments.

When ``array is False``, the value is applied to the fundamental mode
cavities and the frequency of all other cavities is scaled by the same
ratio.

Time lag
~~~~~~~~

The time lag is expressed in values of path lengthening :math:`c\tau`, the 6th
particle coordinate [m].

:py:func:`time_lag = ring.get_rf_timelag(cavpts=None, array=False) <.get_rf_timelag>`

The time lag of the fundamental mode is returned.

:py:func:`ring.set_rf_timelag(time_lag, cavpts=None, array=False, copy=False) <.set_rf_timelag>`

When ``array is False``, the time lag is applied to the fundamental mode
cavities and the time lag of all the other selected cavities is shifted
by the same amount.

All-in-one method
~~~~~~~~~~~~~~~~~

:py:func:`ring.set_cavity(ring, Voltage=None, Frequency=None, TimeLag=None, cavpts=None, array=False, copy=False) <.set_cavity>`

This method sets only the explicitly provided values, the other ones are
left unchanged. For the frequency, a special value :py:class:`at.Frf.NOMINAL <.Frf>`
means nominal frequency, according to the revolution frequency and
harmonic number.

The behaviour of the ``cavpts`` and ``array`` keywords is the same as
for individual methods.

Lattice properties
------------------

The properties provide an even easier way to control the cavities, but
are restricted to the default behaviour of the equivalent Lattice
method:

- cavities are selected by the ``Lattice.cavpts`` attribute if present (lowest frequency
  cavities by default),
- Setting a property modifies the ring in-place (no copy).

:py:attr:`~.Lattice.rf_voltage`

:py:attr:`~.Lattice.rf_frequency`

The special value :py:class:`at.Frf.NOMINAL <.Frf>` means nominal frequency.

:py:attr:`~.Lattice.harmonic_number`

:py:attr:`~.Lattice.rf_timelag`
