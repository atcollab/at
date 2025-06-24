.. _radiation_module:

Radiation
=========

.. rubric:: Functions


.. list-table::

   * - :func:`DipoleRadiation`
     - Compute the radiation integrals in dipoles
   * - :func:`ElementRadiation`
     - - Compute the radiation integrals in dipoles
   * - :func:`ElossRadiation`
     - Compute the radiation integrals in EnergyLoss elements
   * - :func:`FDW`
     - calculate radiation diffusion matrix of a wiggler element.
   * - :func:`WigglerRadiation`
     - Compute the radiation integrals in wigglers
   * - :func:`atdiffmat`
     - quantumDiff    Compute the radiation-diffusion matrix
   * - :func:`atdisable_6d`
     - switches radiation and cavity off
   * - :func:`atenable_6d`
     - switches RF and radiation on
   * - :func:`atenergy`
     - Gets the lattice energy
   * - :func:`atgetU0`
     - Computes Energy loss per turn in eV .
   * - :func:`atradoff`
     - Obsolete: switches RF and radiation off
   * - :func:`atradon`
     - Obsolete: switches RF and radiation on
   * - :func:`atsetenergy`
     - sets the Energy field in all elements.
   * - :func:`attapering`
     - Scale magnet strengths
   * - :func:`check_6d`
     - Checks the presence of longitudinal motion in a lattice.
   * - :func:`check_radiation`
     - Obsolete: check the radiation state of a ring
   * - :func:`findelemraddiffm`
     - 
   * - :func:`findmpoleraddiffmatrix`
     - calculate radiation diffusion matrix of a multipole element
   * - :func:`findthickmpoleraddiffm`
     - 
   * - :func:`findthinmpoleraddiffm`
     - 
   * - :func:`ohmienvelope`
     - calculates equilibrium beam envelope in a
   * - :func:`quantumDiff`
     - Compute the radiation-diffusion matrix
   * - :func:`thickmpoleraddiffm`
     - FIND
   * - :func:`thinmpoleraddiffm`
     - FIND

.. py:function:: DipoleRadiation

   | Compute the radiation integrals in dipoles

.. py:function:: ElementRadiation

   | - Compute the radiation integrals in dipoles
   | and quadrupoles

.. py:function:: ElossRadiation(ring,lindata)

   | Compute the radiation integrals in EnergyLoss elements
   
   | **[i1,i2,i3,i4,i5] = ElossRadiation(ring,lindata)**
   
   | RING       Lattice structure
   | LINDATA    Output of atlinopt for all lattice elements (not used)

.. py:function:: FDW(elem,orbit_in,energy)

   | calculate radiation diffusion matrix of a wiggler element.
   | **diff=FDW(elem,orbit_in,energy)**
   
   | ELEM:      AT wiggler element
   | ORBIT_IN:  input closed orbit
   | ENERGY:    ring energy [GeV]
   
   | **diff=FDW(elem,orbit_in)**
   |    takes energy from the 'Energy' field of the element
   
   |  for use in Ohmi's beam envelope formalism [1]
   |  [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
   | See also :func:`ohmienvelope`

.. py:function:: WigglerRadiation(ring,lindata)

   | Compute the radiation integrals in wigglers
   
   | **[i1,i2,i3,i4,i5] = WigglerRadiation(ring,lindata)**
   
   | RING       Lattice structure
   | LINDATA    Output of atlinopt for all lattice elements
   
   | **WigglerRadiation** computes the radiation integrals for all wigglers with
   | the following approximations:
   
   | - The self-induced dispersion is neglected in I4 and I5, but is is used as
   |   a lower limit for the I5 contribution
   
   |  I1, I2 are integrated analytically
   |  I3 is integrated analytically for a single harmonic, numerically otherwise

.. py:function:: atdiffmat(ring)

   | quantumDiff    Compute the radiation-diffusion matrix
   
   | **[bcum,bs]=atdiffmat(ring)**
   |    RING:       Closed ring AT structure, containing radiative elements and
   |                RF cavity. Radiative elements are identified by a
   |                PassMethod ending with 'RadPass'.
   
   |    BCUM:       Cumulated diffusion matrix
   |    BS:         Cumulative diffusion matrix at the beginning of each element
   
   | **[bcum,bs]=atdiffmat(ring,'orbit',orbitin)**
   |    ORBITIN:    Initial 6-D closed orbit.
   |                In this mode, RING may be a section of a ring.

.. py:function:: atdisable_6d(ring,cavipass,bendpass,quadpass)

   | switches radiation and cavity off
   
   |  **[newring,radindex,cavindex] = atdisable_6d(ring,cavipass,bendpass,quadpass)**
   |     Changes passmethods to turn off cavities, radiation damping and all
   |     elements acting on the particle momentum.
   
   |  The default is to turn everything OFF.,
   
   |   INPUTS:
   |   1. RING      initial AT structure
   |   2. CAVIPASS  pass method for cavities
   |                '' makes no change,
   |                'auto' sets'IdentityPass' or 'DriftPass' depending of cavity length)
   |                anything else is used as the new PassMethod.
   |   3. BENDPASS  pass method for bending magnets
   |                '' makes no change,
   |                'auto' substitutes 'RadPass' with 'Pass' in any method
   |                anything else is used as the new PassMethod.
   |   4. QUADPASS  pass method for quadrupoles
   |                '' makes no change,
   |                'auto' substitutes 'RadPass' with 'Pass' in any method
   |                anything else is used as the new PassMethod.
   
   |   **[...] = atdisable_6d(...[,keyword,value]...)**
   |    The following keywords trigger the processing of the following elements:
   
   |    'allpass'        Defines the default pass method for all elements not
   |                     explicitly specified. Replaces the following default
   |                     values.
   |    'cavipass'       pass method for RF cavities. Default 'auto'
   |    'bendpass'       pass method for bending magnets. Default 'auto'
   |    'quadpass'       pass method for quadrupoles. Default 'auto'
   |    'sextupass'      pass method for sextupoles. Default 'auto'
   |    'octupass'       pass method for bending magnets. Default 'auto'
   |    'multipolepass'  pass method for multipole magnets. Default 'auto'
   |    'wigglerpass'	 pass method for wigglers. Default 'auto'
   |    'quantdiffpass'  pass method for quantum diffusion. Default 'auto'
   |    'energylosspass' pass method for atenergyloss element. Default 'auto'
   |    'simplequantdiffpass' pass method for SimpleQuantDiff element. Default 'auto'
   |    'simpleradiationpass' pass method for SimpleRadiation element. Default 'auto'
   
   |    OUPUTS:
   |    1. NEWRING   Output ring
   |    2. RADINDEX  Indices of elements with radiation
   |    3. CAVINDEX  Indices of active cavities
   
   |   EXAMPLES:
   
   | >> **ringrad=atdisable_6d(ring)**;
   |    Turns off all elements acting on momentum.
   
   | >> **ringrad=atdisable_6d(ring,'auto','allpass','')**;
   |    Turns cavities off and leaves everything else unchanged.
   
   | >> **ringrad=atdisable_6d(ring,'allpass','auto','cavipass','')**;
   |    Turns off everything except RF cavities.
   
   | See also :func:`atenable_6d`, :func:`check_6d`, :func:`atcavityon`, :func:`atcavityoff`

.. py:function:: atenable_6d(ring,cavipass,bendpass,quadpass)

   | switches RF and radiation on
   
   | **[newring,radindex,cavindex] = atenable_6d(ring,cavipass,bendpass,quadpass)**
   |     Changes passmethods to get RF cavity acceleration and radiation
   |     damping.
   
   |  The default is to turn cavities ON and set radiation in dipoles,
   |  quadrupoles and wigglers.
   
   |   INPUTS:
   |   1. RING      initial AT structure
   |   2. CAVIPASS	pass method for cavities
   |                '' makes no change,
   |                'auto' set 'RFCavityPass',
   |                anything else is used as the new PassMethod.
   |   3. BENDPASS	pass method for bending magnets
   |                '' makes no change,
   |                'auto' substitutes 'Pass' with 'RadPass' in any method,
   |                anything else is used as the new PassMethod.
   |   4. QUADPASS	pass method for quadrupoles
   |                '' makes no change,
   |                'auto' substitutes 'Pass' with 'RadPass' in any method,
   |                anything else is used as the new PassMethod.
   
   |   **[...] = atenable_6d(...,keyword,value)**
   |    The following keywords trigger the processing of the following elements:
   
   |    'allpass'        Defines the default pass method for all elements not
   |                     explicitly specified. Replaces the following default
   |                     values.
   |    'cavipass'       pass method for RF cavities. Default 'auto'
   |    'bendpass'       pass method for bending magnets. Default 'auto'
   |    'quadpass'       pass method for quadrupoles. Default 'auto'
   |    'sextupass'      pass method for sextupoles. Default ''
   |    'octupass'       pass method for octupoles. Default ''
   |    'multipolepass'  pass method for multipole magnets. Default ''
   |    'wigglerpass'    pass method for wigglers. Default 'auto'
   |    'quantdiffpass'  pass method for quantum radiation. default 'auto'
   |    'energylosspass' pass method for energyloss element. default 'auto'
   |    'simplequantdiffpass' pass method for SimpleQuantDiff element. Default 'auto'
   |    'simpleradiationpass' pass method for SimpleRadiation element. Default 'auto'
   
   |   OUPUTS:
   |   1. NEWRING   Output ring
   |   2. RADINDEX  Indices of elements with radiation
   |   3. CAVINDEX  Indices of active cavities
   
   |   EXAMPLES:
   
   | >> **ringrad=atenable_6d(ring)**;
   |    Turns cavities on and sets radiation in bending magnets, quadrupoles, energyloss elements, and wigglers (default)
   
   | >> **ringrad=atenable_6d(ring,'auto','allpass','')**;
   |    Turns cavities on and leaves everything else unchanged
   
   | >> **ringrad=atenable_6d(ring,'allpass','','bendpass','auto')**;
   |    Turns on radiation in bending magnets and leaves everything else unchanged
   
   | See also :func:`atdisable_6d`, :func:`check_6d`, :func:`atcavityon`, :func:`atcavityoff`

.. py:function:: atenergy(ring)

   | Gets the lattice energy
   
   |   **energy=atenergy(ring)**
   |   **[energy,periods]=atenergy(ring)**
   |   **[energy,periods,voltage,harmnumber]=atenergy(ring)**
   |   **[energy,periods,voltage,harmnumber,u0]=atenergy(ring)**
   
   |  Warning: To get ENERGY, PERIODS and HARMNUMBER, use atGetRingProperties
   |           To get U0, use atgetU0
   
   |    RING        Ring structure
   
   |    ENERGY      Ring energy
   |        **atenergy** looks for the machine energy in:
   |            1) the 1st 'RingParam' element
   |            2) the 'RFCavity' with the lowest frequency
   |            3) the field "E0" of the global variable "GLOBVAL"
   |            4) The field "Energy" in any element
   |    PERIODS     Number of periods
   |    VOLTAGE     Total RF voltage for the main cavities. The main cavities
   |                are the ones with the lowest frequency
   |    HARMNUMBER  Harmonic number. Computed from the frequency of the main cavities
   |    U0          Total energy loss per turn
   
   | See also :func:`atGetRingProperties`, :func:`atgetU0`, :func:`atsetcavity`

.. py:function:: atgetU0(ring)

   | Computes Energy loss per turn in eV .
   
   | **u0=atgetU0(ring)**   Return the energy loss/turn in eV for the full ring.
   
   |  RING:     Ring structure
   |  U0:       Energy loss per turn in eV
   
   | **u0=atgetU0(...,'periods',periods)** Select the number of periods
   
   |  PERIODS if the number of periods to take into account (Default: full ring)
   
   | **u0=atgetU0(...,'method',method)**	Choose the method
   
   |  METHOD:   'integral': (default) The losses are obtained from
   |                        Losses = Cgamma / 2pi * EGeV^4 * I2
   |                        Takes into account bending magnets and wigglers.
   |            'tracking': The losses are obtained by  tracking without cavities.
   |                        Needs radiation ON, takes into account all radiating elements.
   
   | See also :func:`ringpara`, :func:`atsetcavity`, :func:`atenergy`

.. py:function:: atradoff

   | Obsolete: switches RF and radiation off
   
   |  Kept for compatibility. The function name is misleading, because the
   |  function acts not only on synchrotron radiation, but more generally on
   |  all elements modifying the longitudinal momentum.
   
   |  <a href="matlab:help atdisable_6d">atdisable_6d</a> is an exact copy of this function and should preferably be
   |  used.
   
   | See also :func:`atdisable_6d`, :func:`atenable_6d`, :func:`check_6d`, :func:`atcavityoff`, :func:`atcavityon`

.. py:function:: atradon

   | Obsolete: switches RF and radiation on
   
   |  Kept for compatibility. The function name is misleading, because the
   |  function acts not only on synchrotron radiation, but more generally on
   |  all elements modifying the longitudinal momentum.
   
   |  <a href="matlab:help atenable_6d">atenable_6d</a> is an exact copy of this function and should preferably be
   |  used.
   
   | See also :func:`atenable_6d`, :func:`atdisable_6d`, :func:`check_6d`, :func:`atcavityon`, :func:`atcavityoff`

.. py:function:: atsetenergy(ring,energy)

   | sets the Energy field in all elements.
   |  If no such field exists, it creates it.
   
   |  **newring = atsetenergy(ring,energy)**
   
   |    ring: an AT ring.
   |    Energy: Value to set the Energy field. Units: eV
   |    newring: new AT ring with Energy field set.
   
   |  Example:
   |    Set the energy of the elements in RING to 3 GeV.
   |    **newring = atsetenergy(ring,3e9)**
   
   
   | See also :func:`atenergy`

.. py:function:: attapering(ring)

   | Scale magnet strengths
   
   | **newring=attapering(ring)**   Scales dipole strengths with local energy to
   |    cancel the closed orbit due to synchrotron radiation.
   
   | **newring=attapering(ring,'multipoles', true)**  Scales also the
   |    multipoles to cancel optics errors. The default is true
   
   | **newring=attapering(ring,'niter',niter)** Performs niter iterations (useful
   |    when multipoles are scaled). Default 1
   

.. py:function:: check_6d(ring) returns the radiation state of ring (true/false)

   | Checks the presence of longitudinal motion in a lattice.
   
   | **is_6d = check_6d(ring) returns the radiation state of ring (true/false)**.
   |   Equivalent to IS_6D=<a href="matlab:help atGetRingProperties">atGetRingProperties</a>(RING,'is_6d')
   
   | **is_6d = check_6d(ring, enable)**
   |   Generates an error if IS_6D is different of ENABLE
   
   | **[is_6d, newring] = check_6d(ring,enable,'force')**
   |   The keyword 'force' overrides any error check, and converts
   |   the RING according to the value of ENABLE.
   |   IS_6D contains the status of the RING before conversion.
   
   | **[is_6d, newring] = check_6d(ring,enable,'strict',strict)**
   |   Default, STRICT=true.
   |   If STRICT is true, a difference btw IS_6D and ENABLE produces an error.
   |   If STRICT is false, a difference btw IS_6D and ENABLE produces a warning,
   |   and the RING is converted according to the value of ENABLE.
   |   IS_6D contains the status of the RING before conversion.
   
   | See also :func:`atGetRingProperties`, :func:`atenable_6d`, :func:`atdisable_6d`

.. py:function:: check_radiation

   | Obsolete: check the radiation state of a ring
   
   |  Kept for compatibility> The function name is misleading, because the
   |  function checks not only the presence of synchrotron radiation, but more
   |  generally of all elements modifying the longitudinal momentum.
   
   |  <a href="matlab:help check_6d">check_6d</a> is an exact copy of this function and should preferably be
   |  used.
   
   | See also :func:`check_6d`, :func:`atGetRingProperties`, :func:`atenable_6d`, :func:`atdisable_6d`

.. py:function:: findelemraddiffm

   

.. py:function:: findmpoleraddiffmatrix(elem,orbit_in,energy)

   | calculate radiation diffusion matrix of a multipole element
   | **diff=findmpoleraddiffmatrix(elem,orbit_in,energy)**
   
   | ELEM:      AT multipole element
   | ORBIT_IN:  input closed orbit
   | ENERGY:    ring energy [eV]
   
   | **diff=findmpoleraddiffmatrix(elem,orbit_in)**
   |    takes energy from the 'Energy' field of the element
   
   |  for use in Ohmi's beam envelope formalism [1]
   |  [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
   | See also :func:`ohmienvelope`

.. py:function:: findthickmpoleraddiffm

   

.. py:function:: findthinmpoleraddiffm

   

.. py:function:: ohmienvelope(ring,radelemindex)

   | calculates equilibrium beam envelope in a
   |  circular accelerator using Ohmi's beam envelope formalism [1].
   |  [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
   
   |  **[envelope, rmsdp, rmsbl] = ohmienvelope(ring,radelemindex)**
   |  **[envelope, rmsdp, rmsbl] = ohmienvelope(ring,radelemindex,refpts)**
   
   |  RING    - an AT ring.
   |  RADELEMINDEX - ignored, kept for compatibility
   |  REFPTS  - reference points along the ring. Default: 1
   
   |  ENVELOPE is a structure with fields
   |  Sigma   - [SIGMA(1); SIGMA(2)] - RMS size [m] along
   |            the principal axis of a tilted ellips
   |            Assuming normal distribution exp(-(Z^2)/(2*SIGMA))
   |  Tilt    - Tilt angle of the XY ellipse [rad]
   |            Positive Tilt corresponds to Corkscrew (right)
   |            rotatiom of XY plane around s-axis
   |  R       - 6-by-6 equilibrium envelope matrix R
   
   |  RMSDP   - RMS momentum spread
   |  RMSBL   - RMS bunch length[m]
   
   |  **[envelope, rmsdp, rmsbl, m66, t, orbit] = ohmienvelope(...)**
   |    Returns in addition the 6x6 transfer matrices and the closed orbit
   |    from FINDM66
   
   | See also :func:`atenable_6d`, :func:`findm66`

.. py:function:: quantumDiff(ring)

   | Compute the radiation-diffusion matrix
   
   | **diffmat=quantumDiff(ring)**
   |    RING:       Closed ring AT structure, containing radiative elements and
   |                RF cavity. Radiative elements are identified by a
   |                PassMethod ending with 'RadPass'.
   
   | **diffmat=quantumDiff(line,radindex,orbitin)    (deprecated syntax)**
   | **diffmat=quantumDiff(...,'orbit',orbitin)**
   |    RADINDEX:   Ignored
   |    ORBITIN:    Initial 6-D closed orbit.
   |                In this mode, LINE may be a section of a ring.

.. py:function:: thickmpoleraddiffm

   | FIND

.. py:function:: thinmpoleraddiffm

   | FIND

