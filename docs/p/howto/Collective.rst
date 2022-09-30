Collective
==========

Overview
--------

A collective effects element can be added to any lattice to model
impedance driven collective effects and perform multi-particle tracking.
This element can be constructed with either a user defined table, or by
using the built in functions to construct specific elements (for example
a longitudinal resonator or a transverse resistive wall). The element will
then call the **WakeFieldPass** PassMethod. 

The wake field is applied by first uniformly slicing the full region occupied by the 
particles in the ct co-ordinate. Each particle is then attributed to a
given slice, which is represented by a weight. The mean position in x,y,ct 
is computed for each slice. The kick from one slice to the next (or self kick) can then be computed by taking into account
the differences in offsets of each slice. This kick is then applied to each particle 
within the slice. 

To take into account multi-turn wakes, the wake element has a **TurnHistory** buffer.
Each turn, the mean x,y,ct and weight of each slice is recorded. After each turn, the 
array of ct values is increased by one circumference (to take into account the decay 
between turns). When the kick is computed, the full history of turns is used. 


The package is organised as follows:

**at.collective.wake_functions** contains the analytic wake functions that can be called
by the other classes

The longitudinal resonator wake function is given by [1]

.. math:: W_{z}(\tau) = \left\{ \begin{array}{lr} \alpha R_{s} \;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;   \text{for } \tau=0 \\ 2\alpha R_{s}e^{-\alpha \tau} [\text{cos}(\bar{\omega}\tau) - \frac{\alpha}{\bar{\omega}}\text{sin}(\bar{\omega}\tau)]\;\;\;\; \text{for}\ \tau > 0 \\ \end{array} \right. 

where :math:`R_{s}` is the resonator shunt impedance in Ohms, :math:`\alpha=\omega_{R}/2Q`, :math:`\omega_{R}` is the resonator angular frequency in Hz.rad, :math:`Q` is the resonator quality factor and :math:`\bar{\omega}=\sqrt{\omega_{R}^{2} - \alpha^{2}}`  . The units of the longitudinal wake function are V/C.

The transverse resonator wake function is given by [2]

.. math:: W_{x,y}(z) = -\frac{c R_{s}\omega_{R}}{Q\bar{\omega}}e^{-\alpha z / c} \text{sin}(\bar{\omega} z / c).

but here we use the slightly modified transverse resonator function that is given by

.. math:: W_{x,y}(\tau) = -\frac{ R_{s} \omega_{R}^{2}}{Q\bar{\omega}}e^{-\alpha \tau} \text{sin}(\bar{\omega}\tau).


The definitions are the same as for the longitudinal resonator.

The units of the transverse resonator wake function are V/C/m (but not really)

For both definitions (definitions)

The transverse resistive wall wake function is defined as [3]

.. math:: W_{x,y}(z) = -\frac{2}{\pi b^{3}}\sqrt{\frac{c}{\sigma_{c}}}\frac{1}{|z|^{1/2}}L.

but we use the modified version given by

.. math:: W_{x,y}(\tau) = -\frac{\beta_{x,y}L}{\pi b^{3}}\sqrt{\frac{Z_{0} c}{\pi \sigma_{c} \tau}}.

where :math:`b` is the vacuum chamber half gap in m, :math:`Z_{0}=\pi * 119.9169832` is the impedance of free space in Ohms, c is the velocity of light in m/s, :math:`L` is the length of resistive wall component and :math:`\sigma_{c}` is the conductivity in Siemens. 

The transverse resistive wall function is an approximation, as it clearly diverges for :math:`\tau` close to 0, but when considering multi bunch models the approximation works well. 

Also included is a function to convolve a wake field with a gaussian bunch to compute the wake potential.
This function can be useful to combine multiple wake potentials (for example to include an analytical
wake function with the wake potential output from GDfidL). 

These above functions can be directly called 

.. code:: ipython3

    from at.collective.wake_functions import long_resonator_wf
    from at.collective.wake_functions import transverse_resonator_wf
    from at.collective.wake_functions import transverse_reswall_wf
 
**at.collective.wake_objects** is used to construct a python object which represents a wake field or wake potential. The functions in wake_objects can be called directly to construct and test (and use for purposes outside of tracking in PyAT). The wake object can contain wake fields in multiple planes simultaneously, and may combine wake fields from a table of data with analytical wake functions. All of this is handled by the functions found in this file. 

**at.collective.wake_elements** uses the functions found in **at.collective.wake_objects** to generate a wake field element that can be appended to the AT lattice to be used for tracking. The wake elements can be generated directly, without needing to first make the wake_object.

**at.collective.haissinski** contains functions and methods to solve the Haissinski equation in the presence of a longitudinal wake potential in order to obtain the analytical bunch distribution. 

Generating a Wake Element
-------------------------

We can start with a simple ring. 

.. code:: ipython3

    import at
    ring = at.load_m('at/machine_data/esrf.m')

First we can call the fast_ring function to reduce significantly the number of elements we will need to track

.. code:: ipython3

    fring, _ = at.fast_ring(ring)

First we must define an srange for the wake function. The wake_function will be computed at the values of the srange array, and an interpolation will be made during the tracking if the required dz of the 2 slices falls in between 2 data points. As a way of saving memory, the wake_object contains a useful function for computing the srange such that is is finely sampled only around where the bunches are expected to be. In this example, we will specify how many turns we would like the wake memory to be

.. code:: ipython3

    from at.constants import clight
    from at.collective import Wake

    wturns = 50
    srange_start = 0
    srange_short_end = clight / (2 * ring.get_rf_frequency()) # One half of the bucket width
    sample_fine = 1e-5
    sample_between_bunches = 1e-2   
    bunch_spacing = ring.circumference
    srange_end = wturns * ring.circumference
    
    srange = Wake.build_srange(srange_start, srange_short_end, sample_fine, sample_between_bunches, bunch_spacing, srange_end)
    
Now we can define a longitudinal resonator by calling the LongResonatorElement function from wake_elements. First we need to define some resonator parameters

.. code:: ipython3

    from at.collective.wake_elements import LongResonatorElement

    f_resonator = ring.get_rf_frequency() - 5e4
    qfactor = 4500
    rshunt = 6e6
    current = 0.1   # A
    Nslice = 1
    welem = LongResonatorElement('LongitudinalResonator', ring, srange, f_resonator, qfactor, rshunt, Nturns=wturns, Nslice=Nslice)
    welem.Current = current
    
Finally we can append this to the fast ring

.. code:: ipython3

    fring.append(welem)
    

Using a Wake Table    
------------------

A wake function or wake potential can also be provided from a user defined data or a file. Here we can generate a fake data table using the long_resonator_wf function from at.collective.wake_functions, then we can use it to create a wake element

.. code:: ipython3

    import numpy
    from at.collective import long_resonator_wf
    from at.collective.wake_object import WakeType
    from at.collective.wake_object import WakeComponent
    from at.collective.wake_elements import WakeElement
    
    wf_data = long_resonator_wf(srange, f_resonator, qfactor, rshunt, beta=1)
    
    wa = Wake(srange)
    wa.add(WakeType.TABLE, WakeComponent.Z, srange, wf_data)
    
    welem = WakeElement('wake', ring, wa, Nslice=Nslice)
    
The WakeComponent is used to clearly specify which wake component is being considered. Possible values are Z, DX, DY, QX or QY. 
The WakeType is used to to clearly specify what type of input the add function can expect. Possible values are FILE, RESONATOR, RESWALL or TABLE.
    
Using a Wake File
-----------------

A wake element can also be generated from file. Arguments can be parsed to the add function to describe clearly which columns of the file refer to which parameter. The columns can also be scaled in order to easily sum multiple files or wake contributions.

.. code:: ipython3

    wa = Wake(srange)
    wake_filename = 'filename.txt'

    wa.add(WakeType.FILE, WakeComponent.Z, wake_filename, scol=0, wcol=5, wfact=-1e12)    
    welem = WakeElement('wake', ring, wa, Nslice=Nslice)

Multiple combinations can all be added to one wake element to bring all wake contributions into one wake element

.. code:: ipython3

    wa = Wake(srange)
    wake_filename_z1 = 'filename_z1.txt'
    wf_data_z2 = long_resonator_wf(srange, f_resonator, qfactor, rshunt, beta=1)
    
    wake_filename_dx = 'filename_dx.txt'
    wake_filename_dy = 'filename_dy.txt'

    wa.add(WakeType.FILE, WakeComponent.Z, wake_filename_z1, scol=0, wcol=5, wfact=-1e12)    
    wa.add(WakeType.TABLE, WakeComponent.Z, srange, wf_data_z2)
    wa.add(WakeType.FILE, WakeComponent.DX, wake_filename_dx, scol=0, wcol=1, wfact=1)    
    wa.add(WakeType.FILE, WakeComponent.DY, wake_filename_dy, scol=0, wcol=2, wfact=1)    
    welem = WakeElement('wake', ring, wa, Nslice=Nslice)


Using the Haissinski Class
--------------------------

NOTE: This module is due a re-write and a clean up. But the fundamental process will remain the same.

The Haissinski solver is used to compute the equilibrium beam distribution in the presence of a longitudinal impedance. This class is based entirely on the very nice paper by K. Bane and R. Warnock [4]. In this small overview, we will only talk about how to use it. The details can be seen in the paper of exactly how it is implemented. All the functions within the class are cross referenced with the equations found in the paper. An example file which compares the results of tracking and the results of the Haissinski solver can be found in at/pyat/examples/CollectiveEffects/LongDistribution.py. 

First we initialise a broadband longitudinal resonator wake function in a wake object.

.. code:: ipython3

    from at.collective.wake_object import Wake
    
    circ = 843.977
    freq = 10e9
    qfactor = 1
    Rs = 1e4
    current = 5e-4

    srange = Wake.build_srange(-0.36, 0.36, 1.0e-5, 1.0e-2, circ, circ)

    wobj = Wake.long_resonator(srange, freq, qfactor, rshunt, beta = 1)

Now we need to load and run the Haissinski module. The main parameters here are :math:`m` which defines the number of steps in the distribution, and :math:`k_{max}` which defines the maximum and minimum of the distribution in units of :math:`\sigma_{z}`. numIters is for the number of iterations for the solver to converge to within a convergence criteria of eps. 

.. code:: ipython3

    from at.collective.haissinski import Haissinski

    m = 50 # 30 is quite coarse, 70 or 80 is very fine. 50 is middle
    kmax = 8

    ha = Haissinski(wobj, ring, m=m, kmax=kmax, current=current, numIters = 30, eps=1e-13)
    ha.solve()


The code will now iteratively solve the haissinski equation to determine the beam equilibrium distribution, and will stop running when the distribution no longer changes. Now we can unpack the results and recover some sensible units. 

.. code:: ipython3

    # The x units in the paper are normalised to sigma. So we remove this normalisation.
    ha_x_tmp = ha.q_array*ha.sigma_l 

    # we remove the factor of normalised current
    ha_prof = ha.res/ha.Ic 

    # and now we normalise the profile so the integral is equal to 1
    ha_prof /= numpy.trapz(ha_prof, x=ha_x_tmp) 

    # now we determine the charge center
    ha_cc = numpy.average(ha_x_tmp, weights=ha_prof) 

    # and shift the x position so the bunch is centered around 0
    ha_x = (ha_x_tmp - ha_cc)  

.. image:: haissinski_dist.png



Multi Bunch Collective Effects
------------------------------

Bibliography
------------
[1] A. Chao, 'Physics of Collective Beam Instabilities in High Energy Accelerators', p. 73, Eqn. 2.84

[2] A. Chao, 'Physics of Collective Beam Instabilities in High Energy Accelerators', p. 75, Eqn. 2.88

[3] A. Chao, 'Physics of Collective Beam Instabilities in High Energy Accelerators', p. 59, Eqn. 2.53

[4] "Numerical solution of the Haïssinski equation for the equilibrium state of  a stored electron beam", R. Warnock, K.Bane, Phys. Rev. Acc. and Beams 21, 124401 (2018)



