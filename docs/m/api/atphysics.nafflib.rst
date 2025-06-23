.. _nafflib_module:

nafflib
=======

.. rubric:: Functions


.. list-table::

   * - :func:`calcnaff`
     - Computes NAFF decomposition for a phase space trajectory
   * - :func:`naff_cc`
     - Compile nafflibrary for Matlab
   * - :func:`naff_example`
     - Example to test naff within matlab

.. py:function:: calcnaff(y, yp, win)

   | Computes NAFF decomposition for a phase space trajectory
   |  **[nu amp phase] = calcnaff(y, yp, win)**
   | 
   |   INPUTS
   |   1. Y  - position vector
   |   2. Yp - angle vector
   |   3. WindowType  - Window type - 0 {Default} no windowing
   |                                  1 Window of Hann
   |                                  2 etc
   |   4. nfreq - Maximum number of fundamental frequencies to search for
   |              10 {Default}
   |   5. debug - 1 means debug flag turned on
   |              0 {Default}
   | 
   |   Optional Flags
   |   'Debug' - turn on deubbing flag
   |   'Display' - print ou results
   |   'Hanning' - select Window of Hann, WindowType = 1
   |   'Raw' or 'NoWindow' - select Window of Hann, WindowType = 0
   | 
   |   OUTPUTS
   |   1. frequency - frequency vector with sorted amplitudes
   |                  by default the algorithm try to compute the 10 first fundamental
   |                  frequencies of the system.
   |   2. amplitude - amplitudes associated with fundamental frequencies
   |   3. phase - phases associated with fundamental frequencies
   | 
   |   NOTES
   |   1. Mimimum number of turns is 64 (66)
   |   2. Number of turn has to be a multiple of 6 for internal optimization
   |   reason and just above a power of 2. Example 1026 is divived by 6 and
   |   above 1024 = pow2(10)
   | 
   |   Examples
   |   NT = 9996; % divisible by 6
   |   simple quasiperiodic (even period) motion
   |   y =2+0.1*cos(pi*(0:NT-1))+0.00125*cos(pi/3*(0:NT-1));
   |   yp=2+0.1*sin(pi*(0:NT-1))+0.00125*sin(pi/3*(0:NT-1));
   | 
   |   **[nu ampl phase] = calcnaff(y,yp,1)**; % with windowing

.. py:function:: naff_cc

   | Compile nafflibrary for Matlab
   | 

.. py:function:: naff_example

   | Example to test naff within matlab

