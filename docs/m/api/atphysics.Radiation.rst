.. _radiation_module:

Radiation
=========

.. rubric:: Functions


.. list-table::

   * - :func:`check_6d`
     - CHECK_6D Checks the presence of longitudinal motion in a lattice.
   * - :func:`FDW`
     - FDW calculate radiation diffusion matrix of a wiggler element.
   * - :func:`WigglerRadiation`
     - WIGGLERRADIATION	Compute the radiation integrals in wigglers
   * - :func:`findmpoleraddiffmatrix`
     - FINDMPOLERADDIFFMATRIX calculate radiation diffusion matrix of a multipole element
   * - :func:`atdiffmat`
     - quantumDiff    Compute the radiation-diffusion matrix
   * - :func:`getclass_6d`
     - GETCLASS_6D    Private. Guess class for 6d motion
   * - :func:`atgetU0`
     - ATGETU0 Computes Energy loss per turn in eV .
   * - :func:`findthickmpoleraddiffm`
     - FINDTHICKMPOLERADDIFFM
   * - :func:`check_radiation`
     - CHECK_RADIATION	Obsolete: check the radiation state of a ring
   * - :func:`atenergy`
     - ATENERGY Gets the lattice energy
   * - :func:`radiationoff`
     - RADIATIONOFF turns classical radiation  OFF
   * - :func:`DipoleRadiation`
     - DIPOLERADIATION	Compute the radiation integrals in dipoles
   * - :func:`atenable_6d`
     - ATENABLE_6D switches RF and radiation on
   * - :func:`radiationon`
     - RADIATIONON turns classical radiation  ON
   * - :func:`ElossRadiation`
     - ELOSSRADIATION Compute the radiation integrals in EnergyLoss elements
   * - :func:`ohmienvelope`
     - OHMIENVELOPE calculates equilibrium beam envelope in a
   * - :func:`atradon`
     - ATRADON    Obsolete: switches RF and radiation on
   * - :func:`quantumDiff`
     - QUANTUMDIFF    Compute the radiation-diffusion matrix
   * - :func:`thickmpoleraddiffm`
     - FINDTHICKMPOLERADDIFFM
   * - :func:`attapering`
     - ATTAPERING Scale magnet strengths
   * - :func:`atdisable_6d`
     - ATDISABLE_6D  switches radiation and cavity off
   * - :func:`findelemraddiffm`
     - FINDELEMRADDIFFM
   * - :func:`ElementRadiation`
     - ELEMENTRADIATION - Compute the radiation integrals in dipoles
   * - :func:`atsetenergy`
     - ATSETENERGY sets the Energy field in all elements.
   * - :func:`findthinmpoleraddiffm`
     - FINDTHINMPOLERADDIFFM
   * - :func:`thinmpoleraddiffm`
     - FINDTHINMPOLERADDIFFM
   * - :func:`atradoff`
     - ATRADOFF   Obsolete: switches RF and radiation off

.. py:function:: check_6d

   | 
   | IS_6D = CHECK_6D(RING) Returns the radiation state of RING (true/false).
   |   Equivalent to IS_6D=<a href="matlab:help atGetRingProperties">atGetRingProperties</a>(RING,'is_6d')
   | 
   | IS_6D = CHECK_6D(RING, ENABLE)
   |   Generates an error if IS_6D is different of ENABLE
   | 
   | [IS_6D, NEWRING] = CHECK_6D(RING,ENABLE,'force')
   |   The keyword 'force' overrides any error check, and converts
   |   the RING according to the value of ENABLE.
   |   IS_6D contains the status of the RING before conversion.
   | 
   | [IS_6D, NEWRING] = CHECK_6D(RING,ENABLE,'strict',STRICT)
   |   Default, STRICT=true.
   |   If STRICT is true, a difference btw IS_6D and ENABLE produces an error.
   |   If STRICT is false, a difference btw IS_6D and ENABLE produces a warning,
   |   and the RING is converted according to the value of ENABLE.
   |   IS_6D contains the status of the RING before conversion.
   | 
   |  See also: ATGETRINGPROPERTIES, ATENABLE_6D, ATDISABLE_6D

.. py:function:: FDW

   | DIFF=FDW(ELEM,ORBIT_IN,ENERGY)
   | 
   | ELEM:      AT wiggler element
   | ORBIT_IN:  input closed orbit
   | ENERGY:    ring energy [GeV]
   | 
   | DIFF=FDW(ELEM,ORBIT_IN)
   |    takes energy from the 'Energy' field of the element
   | 
   |  for use in Ohmi's beam envelope formalism [1]
   |  See also OHMIENVELOPE
   |  [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)

.. py:function:: WigglerRadiation

   | 
   | [I1,I2,I3,I4,I5] = WIGGLERRADIATION(RING,LINDATA)
   | 
   | RING       Lattice structure
   | LINDATA    Output of atlinopt for all lattice elements
   | 
   | WigglerRadiation computes the radiation integrals for all wigglers with
   | the following approximations:
   | 
   | - The self-induced dispersion is neglected in I4 and I5, but is is used as
   |   a lower limit for the I5 contribution
   | 
   |  I1, I2 are integrated analytically
   |  I3 is integrated analytically for a single harmonic, numerically otherwise

.. py:function:: findmpoleraddiffmatrix

   | DIFF=FINDMPOLERADDIFFMATRIX(ELEM,ORBIT_IN,ENERGY)
   | 
   | ELEM:      AT multipole element
   | ORBIT_IN:  input closed orbit
   | ENERGY:    ring energy [eV]
   | 
   | DIFF=FINDMPOLERADDIFFMATRIX(ELEM,ORBIT_IN)
   |    takes energy from the 'Energy' field of the element
   | 
   |  for use in Ohmi's beam envelope formalism [1]
   |  See also OHMIENVELOPE
   |  [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)

.. py:function:: atdiffmat

   | 
   | [BCUM,BS]=ATDIFFMAT(RING)
   |    RING:       Closed ring AT structure, containing radiative elements and
   |                RF cavity. Radiative elements are identified by a
   |                PassMethod ending with 'RadPass'.
   | 
   |    BCUM:       Cumulated diffusion matrix
   |    BS:         Cumulative diffusion matrix at the beginning of each element
   | 
   | [BCUM,BS]=ATDIFFMAT(RING,'orbit',ORBITIN)
   |    ORBITIN:    Initial 6-D closed orbit.
   |                In this mode, RING may be a section of a ring.

.. py:function:: getclass_6d

   | 
   |  Returns the element class ignoring its Class field when non-ambiguous
   |  Simplified and faster version of atguessclass

.. py:function:: atgetU0

   | 
   | U0=ATGETU0(RING)   Return the energy loss/turn in eV for the full ring.
   | 
   |  RING:     Ring structure
   |  U0:       Energy loss per turn in eV
   | 
   | U0=ATGETU0(...,'periods',PERIODS) Select the number of periods
   | 
   |  PERIODS if the number of periods to take into account (Default: full ring)
   | 
   | U0=ATGETU0(...,'method',METHOD)	Choose the method
   | 
   |  METHOD:   'integral': (default) The losses are obtained from
   |                        Losses = Cgamma / 2pi * EGeV^4 * I2
   |                        Takes into account bending magnets and wigglers.
   |            'tracking': The losses are obtained by  tracking without cavities.
   |                        Needs radiation ON, takes into account all radiating elements.
   | 
   |  See also RINGPARA ATSETCAVITY ATENERGY

.. py:function:: findthickmpoleraddiffm


.. py:function:: check_radiation

   | 
   |  Kept for compatibility> The function name is misleading, because the
   |  function checks not only the presence of synchrotron radiation, but more
   |  generally of all elements modifying the longitudinal momentum.
   | 
   |  <a href="matlab:help check_6d">check_6d</a> is an exact copy of this function and should preferably be
   |  used.
   | 
   |  See also: CHECK_6D, ATGETRINGPROPERTIES, ATENABLE_6D, ATDISABLE_6D

.. py:function:: atenergy

   | 
   |   ENERGY=ATENERGY(RING)
   |   [ENERGY,PERIODS]=atenergy(RING)
   |   [ENERGY,PERIODS,VOLTAGE,HARMNUMBER]=atenergy(RING)
   |   [ENERGY,PERIODS,VOLTAGE,HARMNUMBER,U0]=atenergy(RING)
   | 
   |  Warning: To get ENERGY, PERIODS and HARMNUMBER, use atGetRingProperties
   |           To get U0, use atgetU0
   | 
   |    RING        Ring structure
   | 
   |    ENERGY      Ring energy
   |        ATENERGY looks for the machine energy in:
   |            1) the 1st 'RingParam' element
   |            2) the 'RFCavity' with the lowest frequency
   |            3) the field "E0" of the global variable "GLOBVAL"
   |            4) The field "Energy" in any element
   |    PERIODS     Number of periods
   |    VOLTAGE     Total RF voltage for the main cavities. The main cavities
   |                are the ones with the lowest frequency
   |    HARMNUMBER  Harmonic number. Computed from the frequency of the main cavities
   |    U0          Total energy loss per turn
   | 
   |   See also atGetRingProperties atgetU0 atsetcavity

.. py:function:: radiationoff

   |   Switch all magnets currently set to use pass-methods
   |   'BndMPoleSymplectic4RadPass' and  'StrMPoleSymplectic4RadPass'
   |   to their equivalents without radiation
   |   'BndMPoleSymplectic4Pass' and  'StrMPoleSymplectic4Pass'
   | 	
   |   NOTES:
   |     1. Deprecated function, use atradoff instead
   | 
   |    See also RADIATIONON, CAVITYON, CAVITYOFF, ATRADON, ATRADOFF

.. py:function:: DipoleRadiation


.. py:function:: atenable_6d

   | 
   | [NEWRING,RADINDEX,CAVINDEX] = ATENABLE_6D(RING,CAVIPASS,BENDPASS,QUADPASS)
   |     Changes passmethods to get RF cavity acceleration and radiation
   |     damping.
   | 
   |  The default is to turn cavities ON and set radiation in dipoles,
   |  quadrupoles and wigglers.
   | 
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
   | 
   |   [...] = ATENABLE_6D(...,keyword,value)
   |    The following keywords trigger the processing of the following elements:
   | 
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
   | 
   |   OUPUTS:
   |   1. NEWRING   Output ring
   |   2. RADINDEX  Indices of elements with radiation
   |   3. CAVINDEX  Indices of active cavities
   | 
   |   EXAMPLES:
   | 
   | >> ringrad=atenable_6d(ring);
   |    Turns cavities on and sets radiation in bending magnets, quadrupoles, energyloss elements, and wigglers (default)
   | 
   | >> ringrad=atenable_6d(ring,'auto','allpass','');
   |    Turns cavities on and leaves everything else unchanged
   | 
   | >> ringrad=atenable_6d(ring,'allpass','','bendpass','auto');
   |    Turns on radiation in bending magnets and leaves everything else unchanged
   | 
   |   See also ATDISABLE_6D, CHECK_6D, ATCAVITYON, ATCAVITYOFF

.. py:function:: radiationon

   | 
   |   Switch all magnets currently set to use pass-methods
   |   'BndMPoleSymplectic4Pass' and  'StrMPoleSymplectic4Pass'
   |   to their equivalents with classical radiation
   |   'BndMPoleSymplectic4RadPass' and  'StrMPoleSymplectic4RadPass'
   | 
   |   NOTES:
   |     1. Deprecated function, use atradon instead
   | 	
   |    See also RADIATIONOFF, CAVITYON, CAVITYOFF, ATRADON, ATRADOFF

.. py:function:: ElossRadiation

   | 
   | [I1,I2,I3,I4,I5] = ELOSSRADIATION(RING,LINDATA)
   | 
   | RING       Lattice structure
   | LINDATA    Output of atlinopt for all lattice elements (not used)

.. py:function:: ohmienvelope

   |  circular accelerator using Ohmi's beam envelope formalism [1].
   |  [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
   | 
   |  [ENVELOPE, RMSDP, RMSBL] = OHMIENVELOPE(RING,RADELEMINDEX)
   |  [ENVELOPE, RMSDP, RMSBL] = OHMIENVELOPE(RING,RADELEMINDEX,REFPTS)
   | 
   |  RING    - an AT ring.
   |  RADELEMINDEX - ignored, kept for compatibility
   |  REFPTS  - reference points along the ring. Default: 1
   | 
   |  ENVELOPE is a structure with fields
   |  Sigma   - [SIGMA(1); SIGMA(2)] - RMS size [m] along
   |            the principal axis of a tilted ellips
   |            Assuming normal distribution exp(-(Z^2)/(2*SIGMA))
   |  Tilt    - Tilt angle of the XY ellipse [rad]
   |            Positive Tilt corresponds to Corkscrew (right)
   |            rotatiom of XY plane around s-axis
   |  R       - 6-by-6 equilibrium envelope matrix R
   | 
   |  RMSDP   - RMS momentum spread
   |  RMSBL   - RMS bunch length[m]
   | 
   |  [ENVELOPE, RMSDP, RMSBL, M66, T, ORBIT] = OHMIENVELOPE(...)
   |    Returns in addition the 6x6 transfer matrices and the closed orbit
   |    from FINDM66
   | 
   |  See also ATENABLE_6D FINDM66

.. py:function:: atradon

   | 
   |  Kept for compatibility. The function name is misleading, because the
   |  function acts not only on synchrotron radiation, but more generally on
   |  all elements modifying the longitudinal momentum.
   | 
   |  <a href="matlab:help atenable_6d">atenable_6d</a> is an exact copy of this function and should preferably be
   |  used.
   | 
   |   See also ATENABLE_6D, ATDISABLE_6D, CHECK_6D, ATCAVITYON, ATCAVITYOFF

.. py:function:: quantumDiff

   | 
   | DIFFMAT=QUANTUMDIFF(RING)
   |    RING:       Closed ring AT structure, containing radiative elements and
   |                RF cavity. Radiative elements are identified by a
   |                PassMethod ending with 'RadPass'.
   | 
   | DIFFMAT=QUANTUMDIFF(LINE,RADINDEX,ORBITIN)    (Deprecated syntax)
   | DIFFMAT=QUANTUMDIFF(...,'orbit',ORBITIN)
   |    RADINDEX:   Ignored
   |    ORBITIN:    Initial 6-D closed orbit.
   |                In this mode, LINE may be a section of a ring.

.. py:function:: thickmpoleraddiffm


.. py:function:: attapering

   | 
   | NEWRING=ATTAPERING(RING)   Scales dipole strengths with local energy to
   |    cancel the closed orbit due to synchrotron radiation.
   | 
   | NEWRING=ATTAPERING(RING,'multipoles', true)  Scales also the
   |    multipoles to cancel optics errors. The default is true
   | 
   | NEWRING=ATTAPERING(RING,'niter',niter) Performs niter iterations (useful
   |    when multipoles are scaled). Default 1
   | 

.. py:function:: atdisable_6d

   | 
   |  [NEWRING,RADINDEX,CAVINDEX] = ATDISABLE_6D(RING,CAVIPASS,BENDPASS,QUADPASS)
   |     Changes passmethods to turn off cavities, radiation damping and all
   |     elements acting on the particle momentum.
   | 
   |  The default is to turn everything OFF., 
   | 
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
   | 
   |   [...] = ATDISABLE_6D(...[,keyword,value]...)
   |    The following keywords trigger the processing of the following elements:
   | 
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
   | 
   |    OUPUTS:
   |    1. NEWRING   Output ring
   |    2. RADINDEX  Indices of elements with radiation
   |    3. CAVINDEX  Indices of active cavities
   | 
   |   EXAMPLES:
   | 
   | >> ringrad=atdisable_6d(ring);
   |    Turns off all elements acting on momentum.
   | 
   | >> ringrad=atdisable_6d(ring,'auto','allpass','');
   |    Turns cavities off and leaves everything else unchanged.
   | 
   | >> ringrad=atdisable_6d(ring,'allpass','auto','cavipass','');
   |    Turns off everything except RF cavities.
   | 
   |   See also ATENABLE_6D, CHECK_6D, ATCAVITYON, ATCAVITYOFF

.. py:function:: findelemraddiffm


.. py:function:: ElementRadiation

   | and quadrupoles

.. py:function:: atsetenergy

   |  If no such field exists, it creates it.
   | 
   |  newring = atsetenergy(ring,Energy)
   | 
   |    ring: an AT ring.
   |    Energy: Value to set the Energy field. Units: eV
   |    newring: new AT ring with Energy field set.
   | 
   |  Example:
   |    Set the energy of the elements in RING to 3 GeV.
   |    NEWRING = atsetenergy(RING,3e9)
   | 
   |  See also atenergy
   | 

.. py:function:: findthinmpoleraddiffm


.. py:function:: thinmpoleraddiffm


.. py:function:: atradoff

   | 
   |  Kept for compatibility. The function name is misleading, because the
   |  function acts not only on synchrotron radiation, but more generally on
   |  all elements modifying the longitudinal momentum.
   | 
   |  <a href="matlab:help atdisable_6d">atdisable_6d</a> is an exact copy of this function and should preferably be
   |  used.
   | 
   |   See also ATDISABLE_6D, ATENABLE_6D, CHECK_6D, ATCAVITYOFF, ATCAVITYON

