.. _element_creation_module:

element_creation
================

.. py:module:: lattice.element_creation

   Create lattice elements

.. rubric:: Functions


.. list-table::

   * - :func:`atM66`
     - Create an element applying an arbitrary 6x6 transfer matrix
   * - :func:`atM66Tijk`
     - ATM66(FAMNAME,M66,Tijk,PASSMETHOD)
   * - :func:`atQuantDiff`
     - creates a quantum diffusion element.
   * - :func:`atSimpleQuantDiff`
     - SimpleQuantDiff creates a simple quantum difusion element
   * - :func:`ataperture`
     - Creates a aperture element
   * - :func:`atbaselem`
     - Create an AT element structure + various checks
   * - :func:`atcorrector`
     - Creates a drift space element with class 'Corrector'
   * - :func:`atcrabcavity`
     - Creates an crab cavity element with Class 'CrabCavity'
   * - :func:`atdampMatElem`
     - creates an element that applies the global damping matrix
   * - :func:`atdrift`
     - Creates a drift space element with Class 'Drift'
   * - :func:`atenergyloss`
     - creates an energy loss element
   * - :func:`atidtable`
     - Creates an ID element
   * - :func:`atinsertiondevicekickmap`
     - creates an insertion device kick-map element
   * - :func:`atmarker`
     - Creates a marker space element
   * - :func:`atmonitor`
     - Creates a Beam Position Monitor element with Class 'Monitor'
   * - :func:`atmultipole`
     - Creates a multipole element
   * - :func:`atquadrupole`
     - Creates a quadrupole element with Class 'Quadrupole'
   * - :func:`atrbend`
     - Creates a rectangular bending magnet element with class 'Bend'
   * - :func:`atrbendtune`
     - Set X0ref and RefDZ for rectangular bending magnets
   * - :func:`atrfcavity`
     - Creates an rfcavity element with Class 'RFCavity'
   * - :func:`atringparam`
     - Creates a RingParameter Element which should go at the beginning of the ring
   * - :func:`atsbend`
     - Creates a sector bending magnet element with class 'Bend'
   * - :func:`atsextupole`
     - Creates a sextupole element with class 'Sextupole'
   * - :func:`atskewquad`
     - Creates a skew quadrupole element with Class 'Multipole'
   * - :func:`atsolenoid`
     - Creates a new solenoid element with Class 'Solenoid'
   * - :func:`atthinmultipole`
     - Creates a thin multipole element
   * - :func:`atvariablemultipole`
     - Creates a variable thin multipole element
   * - :func:`atwiggler`
     - Creates a wiggler

.. py:function:: atM66(famname,atstruct,passmethod)

   | Create an element applying an arbitrary 6x6 transfer matrix
   
   | FAMNAME	family name
   | M66        transfer matrix, defaults to Identity(6)]
   | PASSMETHOD	tracking function, defaults to 'Matrix66Pass'
   
   | **atM66(famname,atstruct,passmethod)**
   |    **atM66** will generate the matrix by calling FINDM66(ATSTRUCT)
   
   | ATSTRUCT   AT structure

.. py:function:: atM66Tijk

   | ATM66(FAMNAME,M66,Tijk,PASSMETHOD)
   |    atM66 creates an element that applies an arbitrary matrix m66
   
   | FAMNAME	family name
   | M66        transfer matrix, defaults to Identity(6)]
   | Tijk       2nd order transfer matrix, defaults to zeros(6,6,6)]
   | PASSMETHOD	tracking function, defaults to 'MatrixTijkPass'
   
   | ATM66(FAMNAME,ATSTRUCT,PASSMETHOD)
   |    atM66 will generate the matrix by calling FINDM66(ATSTRUCT)
   
   | ATSTRUCT   AT structure

.. py:function:: atQuantDiff(famname,diffmat)

   | creates a quantum diffusion element.
   
   | **elem=atQuantDiff(famname,diffmat)** uses the given diffusion matrix
   |    FAMNAME:   family name
   |    DIFFMAT:   Diffusion matrix
   
   | **elem=atQuantDiff(famnane,ring)** computes the diffusion matrix of the ring
   |    FAMNAME:   family name
   |    RING:      lattice with radiation on
   
   | **elem=atQuantDiff(famnane,ring,'orbit0',orbit)** computes the diffusion
   |    matrix of the ring without computing the closed orbit
   |    ORBIT:	closed orbit at beginning of the ring
   |            (this option is useful for the islands)
   
   | The default pass method is 'QuantDiffPass' which uses a global
   | pseudo-random pcg32 stream inplemented in C. More details at:
   | https://github.com/atcollab/at/discussions/879
   
   | See also :func:`quantumDiff`

.. py:function:: atSimpleQuantDiff(famname,...)

   | SimpleQuantDiff creates a simple quantum difusion element
   
   | **elem=atSimpleQuantDiff(famname,...)**
   |    FAMNAME:   family name
   
   | **elem=atSimpleQuantDiff(famname,...,'betax',betax,...)**
   |    BETAX:   Horizontal beta function. Default: 1.0
   
   | **elem=atSimpleQuantDiff(famname,...,'betay',betay,...)**
   |    BETAY:   Vertical beta function. Default: 1.0
   
   | **elem=atSimpleQuantDiff(famname,...,'emitx',emitx,...)**
   |    EMITX:   Horizontal equilibrium emittance. Default: 0.0
   
   | **elem=atSimpleQuantDiff(famname,...,'emity',emity,...)**
   |    EMITY:   Vertical equilibrium emittance. Default: 0.0
   
   | **elem=atSimpleQuantDiff(famname,...,'espread',espread,...)**
   |    ESPREAD: Equilibrium momentum spread. Default: 0.0
   
   | **elem=atSimpleQuantDiff(famname,...,'taux',tau_x,...)**
   |    TAU_X: Horizontal damping time. Default: 0.0
   
   | **elem=atSimpleQuantDiff(famname,...,'tauy',tau_y,...)**
   |    TAU_Y: Vertical damping time. Default: 0.0
   
   | **elem=atSimpleQuantDiff(famname,...,'tauz',tau_z,...)**
   |    TAU_Z: Longitudinal damping time. Default: 0.0

.. py:function:: ataperture

   | Creates a aperture element
   |  **ataperture**(FAMNAME,LIMITS,PASSMETHOD
   |  define physical aperture element (collimator)
   |  lim=[-x,+x,-y,+y];
   
   |  lim={[-x,+x,-y,+y],[-x,+x,-y,+y],[-x,+x,-y,+y],...};
   |  will generate various aperture elements (one for every set of errors)
   
   
   | See also :func:`SetPhysicalAperture`

.. py:function:: atbaselem(famname,method,'fieldname1',value1,...)

   | Create an AT element structure + various checks
   
   | **elem=atbaselem(famname,method,'fieldname1',value1,...)** create AT element
   |    Create an AT element structure and check the consistence of
   |    PolynomA, PolynomB, MaxOrder and NumIntSteps
   
   |   NOTES
   |     1. length of PolynomA and PolynomB are equal (zero padding)
   |     2. MaxOrder is always lenght(PolynomA) - 1

.. py:function:: atcorrector(famname,length,kick,passmethod)

   | Creates a drift space element with class 'Corrector'
   
   |   **atcorrector(famname,length,kick,passmethod)**
   
   |   INPUTS
   |   1. FAMNAME		family name
   |   2. LENGTH		length [m]
   |   3. KICK        [hor. kick, vert. kick] [rad]
   |   4. PASSMETHOD  tracking function, defaults to 'CorrectorPass'
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |   1. Each pair {'FIELDNAME',VALUE} is added to the element
   
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = CorrectorPass
   |                  where req are mandatory field and opt are optional fields
   
   |            atmultipole, atthinmultipole, atmarker
   | See also :func:`atquadrupole`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`

.. py:function:: atcrabcavity(famname,length,voltages,frequency,harmonicnumber)

   | Creates an crab cavity element with Class 'CrabCavity'
   
   |   **atcrabcavity(famname,length,voltages,frequency,harmonicnumber)**
   
   |   INPUTS
   |    1. FAMNAME	    Family name
   |    2. LENGTH		Length [m]
   |    3. VOLTAGES	    Array [Horizontal, Vertical] Peak voltages [V]
   |    4. FREQUENCY	RF frequency [Hz]
   |    5. HARMNUMBER	Harmonic Number
   
   |   OUTPUTS
   |       1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |     **atcrabcavity(famname,...,passmethod,'fieldname1',value1,...)**
   |     Each pair {'FIELDNAME',VALUE} is added to the element
   
   |           atmultipole, atthinmultipole, atmarker, atcorrector
   | See also :func:`atdrift`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`, :func:`atskewquad`, :func:`atrfcavity`

.. py:function:: atdampMatElem(famname,ring,cavipass,bendpass,quadpass)

   | creates an element that applies the global damping matrix
   | **elem=atdampMatElem(famname,ring,cavipass,bendpass,quadpass)**
   
   | FAMNAME:   family name
   | RING:		initial AT structure, without radiation passmethods
   | CAVIPASS:	pass method for cavities (default RFCavityPass)
   | BENDPASS:	pass method for bending magnets. Special values:
   |            '' makes no change,
   |            'auto' wille substitute 'Pass' with 'RadPass' in any method
   |            (default: 'auto')
   | QUADPASS:	pass method for quadrupoles
   |            '' makes no change,
   |            'auto' wille substitute 'Pass' with 'RadPass' in any method
   |            (default: '')
   

.. py:function:: atdrift(famname,length,passmethod)

   | Creates a drift space element with Class 'Drift'
   | **atdrift(famname,length,passmethod)**
   
   |   INPUTS
   |   1. FAMNAME	   - Family name
   |   2. LENGTH	   - Length [m]
   |   3. PASSMETHOD - Tracking function, defaults to 'DriftPass'
   
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |   1. **atdrift(famname,length,passmethod,'fieldname1',value1,...)**
   |     each pair {'fieldname',value} is added to the element
   

.. py:function:: atenergyloss(famname,eloss,passmethod)

   | creates an energy loss element
   
   | **elem=atenergyloss(famname,eloss,passmethod)**
   |    FAMNAME:    family name
   |    ELOSS:      Energy loss [eV]
   |    PASSMETHOD: Tracking methods, defaults to 'IdentityPass'
   
   | the "energy loss" element is taken into account in ATSUMMARY: it adds damping by
   | contributing to the I2 integral, thus reducing the equilibrium emittance.
   | But it does not generate any diffusion. This makes sense only if the losses
   | summarised in the element occur in non-dispersive locations.

.. py:function:: atidtable

   | Creates an ID element
   
   |  FamName	family name
   |  Nslice	number of slices (1 means the wiggler is represented by a
   |            single kick in the center of the device).
   |  filename	name of file with wiggler tracking tables.
   |  Energy    Energy of the machine, needed for scaling
   |  method    tracking function. Defaults to 'IdTablePass'
   
   |  The tracking table is described in
   |  P. Elleaume, "A new approach to the electron beam dynamics in undulators
   |  and wigglers", EPAC92.
   
   |  returns assigned structure with class 'KickMap'

.. py:function:: atinsertiondevicekickmap

   | creates an insertion device kick-map element
   |  Elem = atinsetiondevicekickmap( fname, ...
   |                                  method, ...
   |                                  filename, ...
   |                                  Normalization_energy, ...
   |                                  Nslice, ...
   |                                  length, ...
   |                                  xkick, ...
   |                                  ykick, ...
   |                                  xkick1, ...
   |                                  ykick1, ...
   |                                  xtable, ...
   |                                  ytable ...
   |                                )
   
   |  fname     family name
   |  method    'IdTablePass'
   |  filename  name of the file used to create the element
   |  Normalization_energy    energy to which the field table was scaled
   |  Nslice    number of slices (1 means the wiggler is represented by a
   |            single kick in the center of the device).
   |  length    length of the element
   |  NumX      number of horizontal points
   |  NumY      number of vertical points
   |  xkick     list of x positions
   |  ykick     list of y positions
   |  xkick1    ---
   |  ykick1    ---
   |  xtable    horizontal plane table
   |  ytable    vertical plane table
   
   |  The tracking method is described in
   |  P. Elleaume, "A new approach to the electron beam dynamics in undulators
   |  and wigglers", EPAC92.
   
   |  Returns an element with Class 'InsertionDeviceKickMap'

.. py:function:: atmarker(famname,passmethod)

   | Creates a marker space element
   
   |   **atmarker(famname,passmethod)**
   
   |   INPUTS
   |   1. FAMNAME	 - Family name
   |   2. PASSMETHOD - Tracking function, defaults to 'IdentityPass'
   
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |   1. **atmarker(famname,passmethod,'fieldname1',value1,...)**
   |     each pair {'fieldname',value} is added to the element
   
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = IdentityPass
   |                  where req are mandatory field and opt are optional fields
   
   |           atthinmultipole, atcorrector
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`, :func:`atskewquad`

.. py:function:: atmonitor(famname,'fieldname1',value1,...)

   | Creates a Beam Position Monitor element with Class 'Monitor'
   
   |   INPUTS
   |   1. fname - Family name
   
   |   **atmonitor(famname,'fieldname1',value1,...)**
   |    Each pair {'FIELDNAME',VALUE} is added to the element
   
   |           atmultipole, atthinmultipole, atmarker, atcorrector
   | See also :func:`atdrift`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`

.. py:function:: atmultipole(famname,length,polynoma,polynomb,passmethod)

   | Creates a multipole element
   
   |   **atmultipole(famname,length,polynoma,polynomb,passmethod)**
   
   |   INPUTS
   |   1. FNAME      - Family name
   |   2. LENGTH     - Length[m]
   |   3. POLYNOMA   - Skew [dipole quad sext oct];
   |   4. POLYNOMB   - Normal [dipole quad sext oct];
   |   5. PASSMETHOD - Tracking function. Defaults to 'StrMPoleSymplectic4Pass'
   
   |   OPTIONS (order does not matter)
   |     R1			 -	6 x 6 rotation matrix at the entrance
   | 	 R2        	 -	6 x 6 rotation matrix at the entrance
   | 	 T1			 - 	6 x 1 translation at entrance
   | 	 T2			 -	6 x 1 translation at exit
   | 	 NumIntSteps -   Number of integration steps
   | 	 MaxOrder    -    Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |     1. **atmultipole(famname,length,polynoma,polynomb,passmethod,'fieldname1',value1,...)**
   |    each pair {'fieldname',value} is added to the element
   
   |           atthinmultipole, atmarker, atcorrector
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`, :func:`atskewquad`

.. py:function:: atquadrupole(famname,length,k,passmethod)

   | Creates a quadrupole element with Class 'Quadrupole'
   
   | **atquadrupole(famname,length,k,passmethod)**
   
   |   INPUTS
   |   1. FAMNAME    - Family name
   |   2. LENGTH     - Length [m]
   |   3. K          - Strength [m-2]
   |   4. PASSMETHOD - Tracking function, defaults to 'StrMPoleSymplectic4Pass'
   
   |   OPTIONS (order does not matter)
   |     R1			 -	6 x 6 rotation matrix at the entrance
   | 	 R2        	 -	6 x 6 rotation matrix at the entrance
   | 	 T1			 -	6 x 1 translation at entrance
   | 	 T2			 -	6 x 1 translation at exit
   | 	 NumIntSteps -   Number of integration steps
   | 	 MaxOrder    -   Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = StrMPoleSymplectic4Pass
   |                  where req are mandatory field and opt are optional fields
   |   2. **atquadrupole(famname,length,k,passmethod,'fieldname1',value1,...)**
   |        each pair {'fieldname',value} is added to the element
   
   |   3. Quadrupole fringe field can be activated at element entrance or exit
   |      with option FringeQuadEntrance/FringeQuadExit=0,1,2
   |      Version 0: no fringe field
   |      Version 1: Lee-Whiting formula
   |      Version 2: Lee-Whiting Elegant-like formula where 5 integral need to
   |      be provided
   
   |           atmultipole, atthinmultipole, atmarker, atcorrector, atringparam
   | See also :func:`atdrift`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`, :func:`atskewquad`

.. py:function:: atrbend(famname,length,bendingangle,k,passmethod)

   | Creates a rectangular bending magnet element with class 'Bend'
   
   |   Two calling methods (that can be combined)
   |   **atrbend(famname,length,bendingangle,k,passmethod)**
   
   |   INPUTS
   |   1. FNAME        - Family name
   |   2. LENGTH       - Length of the arc for an on-energy particle
   |                      [m], default to 0
   |   3. BENDINGANGLE - Total bending angle [rad], defaults to 0
   |   4. K			   - Focusing strength, defaults to 0
   |   5. PASSMETHOD   -Tracking function, defaults to 'BndMPoleSymplectic4Pass'
   
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |   1. **atrbend(famname,length,bendingangle,k,passmethod,'fieldname1',value1,...)**
   |     each pair {'fieldname',value} is added to the element
   
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = BndMPoleSymplectic4Pass
   |                  where req are mandatory field and opt are optional fields
   |   2. Model for BndMPoleSymplectic4Pass (Rad) can be selected with extra
   |             fields
   
   |        FringeBendEntrance/FringeBendExit = 0,1,2,3
   |        Version 0 no dipole fringe fields
   |        Version 1 legacy version Brown First Order (K. Brown. A First and Second Order
   |                   Matrix Theory for the Design of Beam Transport Systems and Charged
   |                   Particle Spectrometers. Internal report, SLAC-75, 1982)
   |        Version 2 SOLEIL close to second order of Brown (J. Bengtsson and M. Meddahi.
   |                  Modeling of Beam Dynamics and Comparison with Measurements for
   |                  the Advanced Light Source. London, UK, 1994.)
   |        Version 3 THOMX (Dipole Fringe Field Effects in the ThomX Ring, J. Zhang and
   |                  A. Loulergue, Proceedings of IPAC2013, Shanghai, China)
   
   |        FringeQuadEntrance/FringeQuadExit = 0,1,2
   |        Version 0 no quadrupole fringe fields
   |        Version 1 Lee-Whiting Formula
   |        Version 2 Linear quadrupole fringe field using the 5 integrant a la
   |                  Elegant
   
   |           atmultipole, atthinmultipole, atmarker, atcorrector
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atsextupole`, :func:`atsbend`, :func:`atskewquad`

.. py:function:: atrbendtune(elem)

   | Set X0ref and RefDZ for rectangular bending magnets
   
   | **newelem=atrbendtune(elem)**
   |    Set the X0ref and RefDZ attributes for rectangular bending magnets
   
   | This function must be called after creating a rectangular bending magnet
   | or after setting its polynomA/B attributes. It will set the correct X0ref
   | and RefDZ attributes to get a zero closed orbit for the reference particle.
   
   | Example:
   
   | >> % Identify the rectangular bends
   | >> rbends=atgetcells(ring,...);
   | >> % Set their correct attributes
   | >> ring(rbends)=cellfun(@**atrbendtune**,ring(rbends),'UniformOutput',false);
   
   | Does nothing if the passmethod is not a rectangular bend passmethod

.. py:function:: atrfcavity(famname,length,voltage,frequency,harmonicnumber,energy,passmethod)

   | Creates an rfcavity element with Class 'RFCavity'
   
   |   **atrfcavity(famname,length,voltage,frequency,harmonicnumber,energy,passmethod)**
   
   |   INPUTS
   |    1. FAMNAME	    Family name
   |    2. LENGTH		Length [m]
   |    3. VOLTAGE	    Peak voltage [V]
   |    4. FREQUENCY	RF frequency [Hz]
   |    5. HARMNUMBER	Harmonic Number
   |    6. ENERGY       Energy [eV]
   |    7. PASSMETHOD	Tracking function, defaults to 'RFCavityPass'
   
   |   OUTPUTS
   |       1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |     **atrfcavity(famname,...,passmethod,'fieldname1',value1,...)**
   |     Each pair {'FIELDNAME',VALUE} is added to the element
   
   |   NOTES
   |       1. Fieldname can be called by calling the passmethod
   |          **[req opt] = atrfcavity**
   |                      where req are mandatory field and opt are optional
   |                      fields
   |           atmultipole, atthinmultipole, atmarker, atcorrector
   | See also :func:`atdrift`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`, :func:`atskewquad`

.. py:function:: atringparam(famname,e0,nbperiods)

   | Creates a RingParameter Element which should go at the beginning of the ring
   
   |   **atringparam(famname,e0,nbperiods)**
   
   |   INPUTS
   |   1. FAMNAME	- Family name which may be used as name of Ring
   |   2. E0        - Energy of electrons
   |   3. NBPERIODS - Periodicity of the ring (1 if ring is already expanded)
   
   |   OUTPUTS
   |   1. elem - RingParam class elem
   
   |           atmultipole, atthinmultipole
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`

.. py:function:: atsbend(famname,length,bendingangle,k,passmethod)

   | Creates a sector bending magnet element with class 'Bend'
   
   |   Two calling methods (that can be combined)
   |   **atsbend(famname,length,bendingangle,k,passmethod)**
   
   |   INPUTS
   |   1. FNAME        	Family name
   |   2. LENGTH        Length of the arc for an on-energy particle
   |                      	[m], default to 0
   |   3. BENDINGANGLE	Total bending angle [rad], defaults to 0
   |   4. K				Focusing strength, defaults to 0
   |   5. PASSMETHOD    Tracking function, defaults to 'BndMPoleSymplectic4Pass'
   
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   **atsbend(famname,length,bendingangle,k,passmethod,'fieldname1',value1,...)**
   |   Each pair {'FIELDNAME',VALUE} is added to the element
   
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |          [req opt] = BndMPoleSymplectic4Pass
   |                      where req are mandatory field and opt are optional
   |                      fields
   |   2. Model for BndMPoleSymplectic4Pass (Rad) can be selected with extra
   |             fields
   
   |        FringeBendEntrance/FringeBendExit = 0,1,2,3
   |        Version 0 no dipole fringe fields
   |        Version 1 legacy version Brown First Order (K. Brown. A First and Second Order
   |                   Matrix Theory for the Design of Beam Transport Systems and Charged
   |                   Particle Spectrometers. Internal report, SLAC-75, 1982)
   |        Version 2 SOLEIL close to second order of Brown (J. Bengtsson and M. Meddahi.
   |                  Modeling of Beam Dynamics and Comparison with Measurements for
   |                  the Advanced Light Source. London, UK, 1994.)
   |        Version 3 THOMX (Dipole Fringe Field Effects in the ThomX Ring, J. Zhang and
   |                  A. Loulergue, Proceedings of IPAC2013, Shanghai, China)
   
   |        FringeQuadEntrance/FringeQuadExit = 0,1,2
   |        Version 0 no quadrupole fringe fields
   |        Version 1 Lee-Whiting Formula
   |        Version 2 Linear quadrupole fringe field using the 5 integrant a la
   |                  Elegant
   
   |           atmultipole, atthinmultipole, atmarker, atcorrector
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atsextupole`, :func:`atrbend`

.. py:function:: atsextupole(famname,length,s,passmethod)

   | Creates a sextupole element with class 'Sextupole'
   
   |   **atsextupole(famname,length,s,passmethod)**
   
   |   INPUTS
   | 	 1. FNAME        	family name
   |     2. LENGTH			length [m]
   |     3. S				strength [m-2]
   |     4. PASSMETHOD     tracking function, defaults to 'StrMPoleSymplectic4Pass'
   
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |       1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |     **atsextupole(famname,length,s,passmethod,'fieldname1',value1,...)**
   |     Each pair {'FIELDNAME',VALUE} is added to the element
   
   |             atrbend,atskewquad, atmultipole, atthinmultipole, atmarker,
   |             atcorrector
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atmultipole`, :func:`atsbend`

.. py:function:: atskewquad(famname,length,qs,passmethod)

   | Creates a skew quadrupole element with Class 'Multipole'
   | **atskewquad(famname,length,qs,passmethod)**
   
   |   INPUTS
   |   1. FAMNAME - Family name
   |   2. LENGTH  - Length [m]
   |   3. Qs      - Skew quad strength [m-2]
   
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |   1.  **atskewquad(fname, l, qs, method)**
   
   |           atmultipole, atthinmultipole, atmarker, atcorrector
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`

.. py:function:: atsolenoid

   | Creates a new solenoid element with Class 'Solenoid'
   
   |    Elem =solenoid('FAMILYNAME',Length [m],KS,'METHOD')
   
   |   INPUTS
   | 	1. FamName		  family name
   | 	2. Length	      length[m]
   | 	3. KS             solenoid strength KS [rad/m]
   | 	4. PassMethod     name of the function to use for tracking
   
   
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = BndMPoleSymplectic4Pass
   |                  where req are mandatory field and opt are optional
   |                  fields
   
   |           atthinmultipole, atmarker, atcorrector
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`, :func:`atskewquad`

.. py:function:: atthinmultipole(famname,polynoma,polynomb,passmethod)

   | Creates a thin multipole element
   
   |  **atthinmultipole(famname,polynoma,polynomb,passmethod)**
   
   |   INPUTS
   | 	 1. FNAME           family name
   | 	 2. POLYNOMA        skew [dipole quad sext oct];
   | 	 3. POLYNOMB        normal [dipole quad sext oct];
   | 	 4. PASSMETHOD      tracking function. Defaults to 'ThinMPolePass'
   
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   
   |   OUTPUTS
   |       1. ELEM - Structure with the AT element
   
   |   EXAMPLES
   |     **atthinmultipole(famname,polynoma,polynomb,passmethod,'fieldname1',value1,...)**
   |     Each pair {'FIELDNAME',VALUE} is added to the element
   
   |   NOTES
   |       1. Fieldname can be called by calling the passmethod
   |          [req opt] = BndMPoleSymplectic4Pass
   |                      where req are mandatory field and opt are optional
   |                      fields
   
   |           ATMULTIPOLE, ATMARKER, ATCORRECTOR
   | See also :func:`atdrift`, :func:`atquadrupole`, :func:`atsextupole`, :func:`atsbend`, :func:`atrbend`, :func:`atskewquad`

.. py:function:: atvariablemultipole(famname,mode,passmethod,[key,value]...)

   | Creates a variable thin multipole element
   
   |   **atvariablemultipole(famname,mode,passmethod,[key,value]...)**
   
   |   INPUTS
   |     FNAME          Family name
   |     MODE           Excitation mode: 'SINE', 'WHITENOISE' or 'ARBITRARY'.
   |                    Default: 'SINE'
   |     PASSMETHOD     Tracking function. Default: 'VariableThinMPolePass'
   
   |   OPTIONS (order does not matter)
   |     AMPLITUDEA     Vector or scalar to define the excitation amplitude for
   |                    PolynomA
   |     AMPLITUDEB     Vector or scalar to define the excitation amplitude for
   |                    PolynomA
   |     FREQUENCYA     Frequency of SINE excitation for PolynomA
   |     FREQUENCYB     Frequency of SINE excitation for PolynomB
   |     PHASEA         Phase of SINE excitation for PolynomA
   |     PHASEB         Phase of SINE excitation for PolynomB
   | 	 MAXORDER       Order of the multipole for a scalar amplitude
   |     SEED           Input seed for the random number generator
   |     FUNCA          ARBITRARY excitation turn-by-turn kick list for PolynomA
   |     FUNCB          ARBITRARY excitation turn-by-turn kick list for PolynomB
   |     PERIODIC       If true (default) the user input kick list is repeated
   |     RAMPS          Vector (t0, t1, t2, t3) in turn number to define the ramping of the excitation
   |                    * t<t0: excitation amlpitude is zero
   |                    * t0<t<t1: exciation amplitude is linearly ramped up
   |                    * t1<t<t2: exciation amplitude is constant
   |                    * t2<t<t3: exciation amplitude is linearly ramped down
   |                    * t3<t: exciation amplitude is zero
   
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   
   |   NOTES
   |     1. For all excitation modes at least one amplitude (A or B) is
   |     required. The default excitation is SINE
   |     2. For SINE excitation modes the FREQUENCY corresponding to the input
   |     AMPLITUDE is required
   |     3. For ARBITRARY excitation modes the FUNC corresponding to the input
   |     AMPLITUDE is required
   
   |   EXAMPLES
   
   |  % Create a sinusoidal dipole with amplitude 0.1 mrad and frequency 1 kHz
   |  >> **atvariablemultipole('acm','sine','amplitudeb',1.e-4,'frequencyb',1.e3)**;
   
   |  % Create a white noise dipole excitation of amplitude 0.1 mrad
   |  >> **atvariablemultipole('acm','whitenoise','amplitudeb',1.e-4)**;

.. py:function:: atwiggler(famname, length, lw, bmax, energy, passmethod)

   | Creates a wiggler
   
   | **elem=atwiggler(famname, length, lw, bmax, energy, passmethod)**
   
   |  FAMNAME       family name
   |  LENGTH        total length
   |  LW            Period length
   |  BMAX          Peak magnetic field [T]
   |  ENERGY        Beam energy [eV]
   |  PASSMETHOD    Tracking function. Default 'GWigSymplecticPass'
   
   | **elem=atwiggler(...,'keyword',value...)**
   
   |  Keywords:
   |  Nstep		number of integration steps per period (default 5)
   |  Nmeth		symplectic integration method, 2nd or 4th order: 2 or 4 (default 4)
   |  By		harmonics of the horizontal wiggler. Default [1;1;0;1;1;0]
   |                6xNH matrix, with NH number of harmonics
   |  Bx		harmonics of the vertical wigglers. Default []
   |                6xNV matrix, with NV number of harmonics
   
   | see also: GWigSymplecticPass

