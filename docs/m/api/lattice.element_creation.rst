.. _element_creation_module:

element_creation
================

.. toctree::
   :hidden:

   lattice.element_creation.private

.. rubric:: Modules


.. list-table::

   * - :ref:`private_module`
     - PRIVATE

.. rubric:: Functions


.. list-table::

   * - :func:`atcrabcavity`
     - ATCRABCAVITY Creates an crab cavity element with Class 'CrabCavity'
   * - :func:`atinsertiondevicekickmap`
     - 
   * - :func:`atsbend`
     - ATSBEND Creates a sector bending magnet element with class 'Bend'
   * - :func:`atbaselem`
     - ATBASELEM  Create an AT element structure + various checks
   * - :func:`atrbendtune`
     - ATRBENDTUNE    Set X0ref and RefDZ for rectangular bending magnets
   * - :func:`atrfcavity`
     - ATRFCAVITY Creates an rfcavity element with Class 'RFCavity'
   * - :func:`atwiggler`
     -  ATWIGGLER Creates a wiggler
   * - :func:`atthinmultipole`
     - ATTHINMULTIPOLE Creates a thin multipole element
   * - :func:`atmonitor`
     - ATMONITOR Creates a Beam Position Monitor element with Class 'Monitor'
   * - :func:`atenergyloss`
     - atenergyloss creates an energy loss element
   * - :func:`ataperture`
     - ATAPERTURE Creates a aperture element
   * - :func:`atsextupole`
     - ATSEXTUPOLE Creates a sextupole element with class 'Sextupole'
   * - :func:`atsolenoid`
     - ATSOLENOID Creates a new solenoid element with Class 'Solenoid'
   * - :func:`atrbend`
     - ATRBEND Creates a rectangular bending magnet element with class 'Bend'
   * - :func:`atquadrupole`
     - ATQUADRUPOLE Creates a quadrupole element with Class 'Quadrupole'
   * - :func:`atM66Tijk`
     - ATM66(FAMNAME,M66,Tijk,PASSMETHOD)
   * - :func:`atidtable`
     -  ATIDTABLE Creates an ID element
   * - :func:`atSimpleQuantDiff`
     - SimpleQuantDiff creates a simple quantum difusion element
   * - :func:`atdampMatElem`
     -    atdampMatElem creates an element that applies the global damping matrix
   * - :func:`atQuantDiff`
     - atQuantDiff creates a quantum diffusion element.
   * - :func:`atmarker`
     - ATMARKER Creates a marker space element
   * - :func:`atcorrector`
     - ATCORRECTOR Creates a drift space element with class 'Corrector'
   * - :func:`atdrift`
     - ATDRIFT Creates a drift space element with Class 'Drift'
   * - :func:`atvariablemultipole`
     - ATVARIABLEMULTIPOLE Creates a variable thin multipole element
   * - :func:`atringparam`
     - ATRINGPARAM Creates a RingParameter Element which should go at the beginning of the ring
   * - :func:`atM66`
     - ATM66  Create an element applying an arbitrary 6x6 transfer matrix
   * - :func:`atmultipole`
     - ATMULTIPOLE Creates a multipole element
   * - :func:`atskewquad`
     - ATSKEWQUAD Creates a skew quadrupole element with Class 'Multipole'

.. py:function:: atcrabcavity

   | ATCRABCAVITY Creates an crab cavity element with Class 'CrabCavity'
   | 
   |   ATCRABCAVITY(FAMNAME,LENGTH,VOLTAGES,FREQUENCY,HARMONICNUMBER)
   | 
   |   INPUTS
   |    1. FAMNAME	    Family name
   |    2. LENGTH		Length [m]
   |    3. VOLTAGES	    Array [Horizontal, Vertical] Peak voltages [V]
   |    4. FREQUENCY	RF frequency [Hz]
   |    5. HARMNUMBER	Harmonic Number
   | 
   |   OUTPUTS
   |       1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |     ATCRABCAVITY(FAMNAME,...,PASSMETHOD,'FIELDNAME1',VALUE1,...)
   |     Each pair {'FIELDNAME',VALUE} is added to the element
   | 
   | See also  atdrift, atsextupole, atsbend, atrbend, atskewquad, atrfcavity
   |           atmultipole, atthinmultipole, atmarker, atcorrector

.. py:function:: atinsertiondevicekickmap


.. py:function:: atsbend

   | ATSBEND Creates a sector bending magnet element with class 'Bend'
   | 
   |   Two calling methods (that can be combined)
   |   ATSBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD)
   | 
   |   INPUTS
   |   1. FNAME        	Family name
   |   2. LENGTH        Length of the arc for an on-energy particle
   |                      	[m], default to 0
   |   3. BENDINGANGLE	Total bending angle [rad], defaults to 0
   |   4. K				Focusing strength, defaults to 0
   |   5. PASSMETHOD    Tracking function, defaults to 'BndMPoleSymplectic4Pass'
   | 
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   ATSBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD,'FIELDNAME1',VALUE1,...)
   |   Each pair {'FIELDNAME',VALUE} is added to the element
   | 
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |          [req opt] = BndMPoleSymplectic4Pass
   |                      where req are mandatory field and opt are optional
   |                      fields
   |   2. Model for BndMPoleSymplectic4Pass (Rad) can be selected with extra
   |             fields
   | 
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
   | 
   |        FringeQuadEntrance/FringeQuadExit = 0,1,2
   |        Version 0 no quadrupole fringe fields
   |        Version 1 Lee-Whiting Formula
   |        Version 2 Linear quadrupole fringe field using the 5 integrant a la
   |                  Elegant
   | 
   |   See also atdrift, atquadrupole, atsextupole, atrbend
   |           atmultipole, atthinmultipole, atmarker, atcorrector

.. py:function:: atbaselem

   | ATBASELEM  Create an AT element structure + various checks
   | 
   | ELEM=ATBASELEM(FAMNAME,METHOD,'FIELDNAME1',VALUE1,...) create AT element
   |    Create an AT element structure and check the consistence of
   |    PolynomA, PolynomB, MaxOrder and NumIntSteps
   | 
   |   NOTES
   |     1. length of PolynomA and PolynomB are equal (zero padding)
   |     2. MaxOrder is always lenght(PolynomA) - 1

.. py:function:: atrbendtune

   | ATRBENDTUNE    Set X0ref and RefDZ for rectangular bending magnets
   | 
   | NEWELEM=ATRBENDTUNE(ELEM)
   |    Set the X0ref and RefDZ attributes for rectangular bending magnets
   | 
   | This function must be called after creating a rectangular bending magnet
   | or after setting its polynomA/B attributes. It will set the correct X0ref
   | and RefDZ attributes to get a zero closed orbit for the reference particle.
   | 
   | Example:
   | 
   | >> % Identify the rectangular bends
   | >> rbends=atgetcells(ring,...);
   | >> % Set their correct attributes
   | >> ring(rbends)=cellfun(@atrbendtune,ring(rbends),'UniformOutput',false);
   | 
   | Does nothing if the passmethod is not a rectangular bend passmethod

.. py:function:: atrfcavity

   | ATRFCAVITY Creates an rfcavity element with Class 'RFCavity'
   | 
   |   ATRFCAVITY(FAMNAME,LENGTH,VOLTAGE,FREQUENCY,HARMONICNUMBER,ENERGY,PASSMETHOD)
   | 
   |   INPUTS
   |    1. FAMNAME	    Family name
   |    2. LENGTH		Length [m]
   |    3. VOLTAGE	    Peak voltage [V]
   |    4. FREQUENCY	RF frequency [Hz]
   |    5. HARMNUMBER	Harmonic Number
   |    6. ENERGY       Energy [eV]
   |    7. PASSMETHOD	Tracking function, defaults to 'RFCavityPass'
   | 
   |   OUTPUTS
   |       1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |     ATRFCAVITY(FAMNAME,...,PASSMETHOD,'FIELDNAME1',VALUE1,...)
   |     Each pair {'FIELDNAME',VALUE} is added to the element
   | 
   |   NOTES
   |       1. Fieldname can be called by calling the passmethod
   |          [req opt] = atrfcavity
   |                      where req are mandatory field and opt are optional
   |                      fields
   | %See also  atdrift, atsextupole, atsbend, atrbend, atskewquad
   |           atmultipole, atthinmultipole, atmarker, atcorrector

.. py:function:: atwiggler

   |  ATWIGGLER Creates a wiggler
   | 
   | ELEM=ATWIGGLER(FAMNAME, LENGTH, LW, BMAX, ENERGY, PASSMETHOD)
   | 
   |  FAMNAME       family name
   |  LENGTH        total length
   |  LW            Period length
   |  BMAX          Peak magnetic field [T]
   |  ENERGY        Beam energy [eV]
   |  PASSMETHOD    Tracking function. Default 'GWigSymplecticPass'
   | 
   | ELEM=ATWIGGLER(...,'keyword',value...)
   | 
   |  Keywords:
   |  Nstep		number of integration steps per period (default 5)
   |  Nmeth		symplectic integration method, 2nd or 4th order: 2 or 4 (default 4)
   |  By		harmonics of the horizontal wiggler. Default [1;1;0;1;1;0]
   |                6xNH matrix, with NH number of harmonics
   |  Bx		harmonics of the vertical wigglers. Default []
   |                6xNV matrix, with NV number of harmonics
   | 
   | see also: GWigSymplecticPass

.. py:function:: atthinmultipole

   | ATTHINMULTIPOLE Creates a thin multipole element
   | 
   |  ATTHINMULTIPOLE(FAMNAME,POLYNOMA,POLYNOMB,PASSMETHOD)
   | 
   |   INPUTS
   | 	 1. FNAME           family name
   | 	 2. POLYNOMA        skew [dipole quad sext oct];
   | 	 3. POLYNOMB        normal [dipole quad sext oct];
   | 	 4. PASSMETHOD      tracking function. Defaults to 'ThinMPolePass'
   | 
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |       1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |     ATTHINMULTIPOLE(FAMNAME,POLYNOMA,POLYNOMB,PASSMETHOD,'FIELDNAME1',VALUE1,...)
   |     Each pair {'FIELDNAME',VALUE} is added to the element
   | 
   |   NOTES
   |       1. Fieldname can be called by calling the passmethod
   |          [req opt] = BndMPoleSymplectic4Pass
   |                      where req are mandatory field and opt are optional
   |                      fields
   | 
   | See also  ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND ATSKEWQUAD,
   |           ATMULTIPOLE, ATMARKER, ATCORRECTOR

.. py:function:: atmonitor

   | ATMONITOR Creates a Beam Position Monitor element with Class 'Monitor'
   | 
   |   INPUTS
   |   1. fname - Family name
   | 
   |   ATMONITOR(FAMNAME,'FIELDNAME1',VALUE1,...)
   |    Each pair {'FIELDNAME',VALUE} is added to the element
   | 
   |   See also atdrift, atsextupole, atsbend, atrbend
   |           atmultipole, atthinmultipole, atmarker, atcorrector

.. py:function:: atenergyloss

   | atenergyloss creates an energy loss element
   | 
   | ELEM=ATENERGYLOSS(FAMNAME,ELOSS,PASSMETHOD)
   |    FAMNAME:    family name
   |    ELOSS:      Energy loss [eV]
   |    PASSMETHOD: Tracking methods, defaults to 'IdentityPass'
   | 
   | the "energy loss" element is taken into account in ATSUMMARY: it adds damping by
   | contributing to the I2 integral, thus reducing the equilibrium emittance.
   | But it does not generate any diffusion. This makes sense only if the losses
   | summarised in the element occur in non-dispersive locations.

.. py:function:: ataperture

   | ATAPERTURE Creates a aperture element
   |  ATAPERTURE(FAMNAME,LIMITS,PASSMETHOD
   |  define physical aperture element (collimator)
   |  lim=[-x,+x,-y,+y];
   | 
   |  lim={[-x,+x,-y,+y],[-x,+x,-y,+y],[-x,+x,-y,+y],...};
   |  will generate various aperture elements (one for every set of errors)
   | 
   |   See also setphysicalaperture
   | 

.. py:function:: atsextupole

   | ATSEXTUPOLE Creates a sextupole element with class 'Sextupole'
   | 
   |   ATSEXTUPOLE(FAMNAME,LENGTH,S,PASSMETHOD)
   | 
   |   INPUTS
   | 	 1. FNAME        	family name
   |     2. LENGTH			length [m]
   |     3. S				strength [m-2]
   |     4. PASSMETHOD     tracking function, defaults to 'StrMPoleSymplectic4Pass'
   | 
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |       1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |     ATSEXTUPOLE(FAMNAME,LENGTH,S,PASSMETHOD,'FIELDNAME1',VALUE1,...)
   |     Each pair {'FIELDNAME',VALUE} is added to the element
   | 
   |   See also: atdrift, atquadrupole, atmultipole, atsbend,
   |             atrbend,atskewquad, atmultipole, atthinmultipole, atmarker,
   |             atcorrector

.. py:function:: atsolenoid

   | ATSOLENOID Creates a new solenoid element with Class 'Solenoid'
   | 
   |    Elem =solenoid('FAMILYNAME',Length [m],KS,'METHOD')
   | 
   |   INPUTS
   | 	1. FamName		  family name
   | 	2. Length	      length[m]
   | 	3. KS             solenoid strength KS [rad/m]
   | 	4. PassMethod     name of the function to use for tracking
   | 
   | 
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = BndMPoleSymplectic4Pass
   |                  where req are mandatory field and opt are optional
   |                  fields
   | 
   |   See also atdrift, atquadrupole, atsextupole, atsbend, atrbend atskewquad,
   |           atthinmultipole, atmarker, atcorrector

.. py:function:: atrbend

   | ATRBEND Creates a rectangular bending magnet element with class 'Bend'
   | 
   |   Two calling methods (that can be combined)
   |   ATRBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD)
   | 
   |   INPUTS
   |   1. FNAME        - Family name
   |   2. LENGTH       - Length of the arc for an on-energy particle
   |                      [m], default to 0
   |   3. BENDINGANGLE - Total bending angle [rad], defaults to 0
   |   4. K			   - Focusing strength, defaults to 0
   |   5. PASSMETHOD   -Tracking function, defaults to 'BndMPoleSymplectic4Pass'
   | 
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |   1. atrbend(famname,length,bendingangle,k,passmethod,'fieldname1',value1,...)
   |     each pair {'fieldname',value} is added to the element
   | 
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = BndMPoleSymplectic4Pass
   |                  where req are mandatory field and opt are optional fields
   |   2. Model for BndMPoleSymplectic4Pass (Rad) can be selected with extra
   |             fields
   | 
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
   | 
   |        FringeQuadEntrance/FringeQuadExit = 0,1,2
   |        Version 0 no quadrupole fringe fields
   |        Version 1 Lee-Whiting Formula
   |        Version 2 Linear quadrupole fringe field using the 5 integrant a la
   |                  Elegant
   | 
   |   See also atdrift, atquadrupole, atsextupole, atsbend, atskewquad,
   |           atmultipole, atthinmultipole, atmarker, atcorrector

.. py:function:: atquadrupole

   | ATQUADRUPOLE Creates a quadrupole element with Class 'Quadrupole'
   | 
   | ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD)
   | 
   |   INPUTS
   |   1. FAMNAME    - Family name
   |   2. LENGTH     - Length [m]
   |   3. K          - Strength [m-2]
   |   4. PASSMETHOD - Tracking function, defaults to 'StrMPoleSymplectic4Pass'
   | 
   |   OPTIONS (order does not matter)
   |     R1			 -	6 x 6 rotation matrix at the entrance
   | 	 R2        	 -	6 x 6 rotation matrix at the entrance
   | 	 T1			 -	6 x 1 translation at entrance
   | 	 T2			 -	6 x 1 translation at exit
   | 	 NumIntSteps -   Number of integration steps
   | 	 MaxOrder    -   Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = StrMPoleSymplectic4Pass
   |                  where req are mandatory field and opt are optional fields
   |   2. atquadrupole(famname,length,k,passmethod,'fieldname1',value1,...)
   |        each pair {'fieldname',value} is added to the element
   | 
   |   3. Quadrupole fringe field can be activated at element entrance or exit
   |      with option FringeQuadEntrance/FringeQuadExit=0,1,2
   |      Version 0: no fringe field
   |      Version 1: Lee-Whiting formula
   |      Version 2: Lee-Whiting Elegant-like formula where 5 integral need to
   |      be provided
   | 
   |   See also atdrift, atsextupole, atsbend, atrbend, atskewquad,
   |           atmultipole, atthinmultipole, atmarker, atcorrector, atringparam

.. py:function:: atM66Tijk

   | ATM66(FAMNAME,M66,Tijk,PASSMETHOD)
   |    atM66 creates an element that applies an arbitrary matrix m66
   | 
   | FAMNAME	family name
   | M66        transfer matrix, defaults to Identity(6)]
   | Tijk       2nd order transfer matrix, defaults to zeros(6,6,6)]
   | PASSMETHOD	tracking function, defaults to 'MatrixTijkPass'
   | 
   | ATM66(FAMNAME,ATSTRUCT,PASSMETHOD)
   |    atM66 will generate the matrix by calling FINDM66(ATSTRUCT)
   | 
   | ATSTRUCT   AT structure

.. py:function:: atidtable

   |  ATIDTABLE Creates an ID element
   | 
   |  FamName	family name
   |  Nslice	number of slices (1 means the wiggler is represented by a
   |            single kick in the center of the device).
   |  filename	name of file with wiggler tracking tables.
   |  Energy    Energy of the machine, needed for scaling
   |  method    tracking function. Defaults to 'IdTablePass'
   | 
   |  The tracking table is described in
   |  P. Elleaume, "A new approach to the electron beam dynamics in undulators
   |  and wigglers", EPAC92.
   | 
   |  returns assigned structure with class 'KickMap'

.. py:function:: atSimpleQuantDiff

   | SimpleQuantDiff creates a simple quantum difusion element
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...)
   |    FAMNAME:   family name
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'betax',BETAX,...)
   |    BETAX:   Horizontal beta function. Default: 1.0
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'betay',BETAY,...)
   |    BETAY:   Vertical beta function. Default: 1.0
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'emitx',EMITX,...)
   |    EMITX:   Horizontal equilibrium emittance. Default: 0.0
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'emity',EMITY,...)
   |    EMITY:   Vertical equilibrium emittance. Default: 0.0
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'espread',ESPREAD,...)
   |    ESPREAD: Equilibrium momentum spread. Default: 0.0
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'taux',TAU_X,...)
   |    TAU_X: Horizontal damping time. Default: 0.0
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'tauy',TAU_Y,...)
   |    TAU_Y: Vertical damping time. Default: 0.0
   | 
   | ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'tauz',TAU_Z,...)
   |    TAU_Z: Longitudinal damping time. Default: 0.0

.. py:function:: atdampMatElem

   |    atdampMatElem creates an element that applies the global damping matrix
   | ELEM=ATDAMPMATELEM(FAMNAME,RING,CAVIPASS,BENDPASS,QUADPASS)
   | 
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
   | 

.. py:function:: atQuantDiff

   | atQuantDiff creates a quantum diffusion element.
   | 
   | ELEM=ATQUANTDIFF(FAMNAME,DIFFMAT) uses the given diffusion matrix
   |    FAMNAME:   family name
   |    DIFFMAT:   Diffusion matrix
   | 
   | ELEM=ATQUANTDIFF(FAMNANE,RING) computes the diffusion matrix of the ring
   |    FAMNAME:   family name
   |    RING:      lattice with radiation on
   | 
   | ELEM=ATQUANTDIFF(FAMNANE,RING,'orbit0',ORBIT) computes the diffusion
   |    matrix of the ring without computing the closed orbit
   |    ORBIT:	closed orbit at beginning of the ring
   |            (this option is useful for the islands)
   | 
   | The default pass method is 'QuantDiffPass' which uses a global
   | pseudo-random pcg32 stream inplemented in C. More details at:
   | https://github.com/atcollab/at/discussions/879
   | 
   | See also quantumDiff

.. py:function:: atmarker

   | ATMARKER Creates a marker space element
   | 
   |   ATMARKER(FAMNAME,PASSMETHOD)
   | 
   |   INPUTS
   |   1. FAMNAME	 - Family name
   |   2. PASSMETHOD - Tracking function, defaults to 'IdentityPass'
   | 
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |   1. atmarker(famname,passmethod,'fieldname1',value1,...)
   |     each pair {'fieldname',value} is added to the element
   | 
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = IdentityPass
   |                  where req are mandatory field and opt are optional fields
   | 
   | See also  atdrift, atquadrupole, atsextupole, atsbend, atrbend atskewquad,
   |           atthinmultipole, atcorrector

.. py:function:: atcorrector

   | ATCORRECTOR Creates a drift space element with class 'Corrector'
   | 
   |   atcorrector(FAMNAME,LENGTH,KICK,PASSMETHOD)
   | 
   |   INPUTS
   |   1. FAMNAME		family name
   |   2. LENGTH		length [m]
   |   3. KICK        [hor. kick, vert. kick] [rad]
   |   4. PASSMETHOD  tracking function, defaults to 'CorrectorPass'
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |   1. Each pair {'FIELDNAME',VALUE} is added to the element
   | 
   |   NOTES
   |   1. Fieldname can be called by calling the passmethod
   |      [req opt] = CorrectorPass
   |                  where req are mandatory field and opt are optional fields
   | 
   |   See also atquadrupole, atsextupole, atsbend, atrbend
   |            atmultipole, atthinmultipole, atmarker

.. py:function:: atdrift

   | ATDRIFT Creates a drift space element with Class 'Drift'
   | ATDRIFT(FAMNAME,LENGTH,PASSMETHOD)
   | 
   |   INPUTS
   |   1. FAMNAME	   - Family name
   |   2. LENGTH	   - Length [m]
   |   3. PASSMETHOD - Tracking function, defaults to 'DriftPass'
   | 
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |   1. atdrift(famname,length,passmethod,'fieldname1',value1,...)
   |     each pair {'fieldname',value} is added to the element
   | 

.. py:function:: atvariablemultipole

   | ATVARIABLEMULTIPOLE Creates a variable thin multipole element
   | 
   |   ATVARIABLEMULTIPOLE(FAMNAME,MODE,PASSMETHOD,[KEY,VALUE]...)
   | 
   |   INPUTS
   |     FNAME          Family name
   |     MODE           Excitation mode: 'SINE', 'WHITENOISE' or 'ARBITRARY'.
   |                    Default: 'SINE'
   |     PASSMETHOD     Tracking function. Default: 'VariableThinMPolePass'
   | 
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
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   NOTES
   |     1. For all excitation modes at least one amplitude (A or B) is
   |     required. The default excitation is SINE
   |     2. For SINE excitation modes the FREQUENCY corresponding to the input
   |     AMPLITUDE is required
   |     3. For ARBITRARY excitation modes the FUNC corresponding to the input
   |     AMPLITUDE is required
   | 
   |   EXAMPLES
   | 
   |  % Create a sinusoidal dipole with amplitude 0.1 mrad and frequency 1 kHz
   |  >> atvariablemultipole('ACM','SINE','AmplitudeB',1.e-4,'FrequencyB',1.e3);
   | 
   |  % Create a white noise dipole excitation of amplitude 0.1 mrad
   |  >> atvariablemultipole('ACM','WHITENOISE','AmplitudeB',1.e-4);

.. py:function:: atringparam

   | ATRINGPARAM Creates a RingParameter Element which should go at the beginning of the ring
   | 
   |   atringparam(FAMNAME,E0,NBPERIODS)
   | 
   |   INPUTS
   |   1. FAMNAME	- Family name which may be used as name of Ring
   |   2. E0        - Energy of electrons
   |   3. NBPERIODS - Periodicity of the ring (1 if ring is already expanded)
   | 
   |   OUTPUTS
   |   1. elem - RingParam class elem
   | 
   |   See also atdrift, atquadrupole, atsextupole, atsbend, atrbend
   |           atmultipole, atthinmultipole

.. py:function:: atM66

   | ATM66  Create an element applying an arbitrary 6x6 transfer matrix
   | 
   | FAMNAME	family name
   | M66        transfer matrix, defaults to Identity(6)]
   | PASSMETHOD	tracking function, defaults to 'Matrix66Pass'
   | 
   | ATM66(FAMNAME,ATSTRUCT,PASSMETHOD)
   |    atM66 will generate the matrix by calling FINDM66(ATSTRUCT)
   | 
   | ATSTRUCT   AT structure

.. py:function:: atmultipole

   | ATMULTIPOLE Creates a multipole element
   | 
   |   ATMULTIPOLE(FAMNAME,LENGTH,POLYNOMA,POLYNOMB,PASSMETHOD)
   | 
   |   INPUTS
   |   1. FNAME      - Family name
   |   2. LENGTH     - Length[m]
   |   3. POLYNOMA   - Skew [dipole quad sext oct];
   |   4. POLYNOMB   - Normal [dipole quad sext oct];
   |   5. PASSMETHOD - Tracking function. Defaults to 'StrMPoleSymplectic4Pass'
   | 
   |   OPTIONS (order does not matter)
   |     R1			 -	6 x 6 rotation matrix at the entrance
   | 	 R2        	 -	6 x 6 rotation matrix at the entrance
   | 	 T1			 - 	6 x 1 translation at entrance
   | 	 T2			 -	6 x 1 translation at exit
   | 	 NumIntSteps -   Number of integration steps
   | 	 MaxOrder    -    Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |     1. atmultipole(famname,length,polynoma,polynomb,passmethod,'fieldname1',value1,...)
   |    each pair {'fieldname',value} is added to the element
   | 
   |   See also atdrift, atquadrupole, atsextupole, atsbend, atrbend, atskewquad,
   |           atthinmultipole, atmarker, atcorrector

.. py:function:: atskewquad

   | ATSKEWQUAD Creates a skew quadrupole element with Class 'Multipole'
   | atskewquad(famname,length,qs,passmethod)
   | 
   |   INPUTS
   |   1. FAMNAME - Family name
   |   2. LENGTH  - Length [m]
   |   3. Qs      - Skew quad strength [m-2]
   | 
   |   OPTIONS (order does not matter)
   |     R1				6 x 6 rotation matrix at the entrance
   | 	 R2        		6 x 6 rotation matrix at the entrance
   | 	 T1				6 x 1 translation at entrance
   | 	 T2				6 x 1 translation at exit
   | 	 NumIntSteps    Number of integration steps
   | 	 MaxOrder       Max Order for multipole (1 up to quadrupole)
   | 
   |   OUTPUTS
   |   1. ELEM - Structure with the AT element
   | 
   |   EXAMPLES
   |   1.  atskewquad(Fname, L, Qs, method)
   | 
   |   See also atdrift, atquadrupole, atsextupole, atsbend, atrbend,
   |           atmultipole, atthinmultipole, atmarker, atcorrector

