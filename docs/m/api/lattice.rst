.. _lattice_module:

lattice
=======

.. toctree::
   :hidden:

   lattice.Converters
   lattice.Paramgroup
   lattice.element_creation
   lattice.survey

.. rubric:: Modules


.. list-table::

   * - :ref:`converters_module`
     - % CONVERTERS

   * - :ref:`paramgroup_module`
     - PARAMGROUP
   * - :ref:`element_creation_module`
     - % ELEMENT_CREATION

   * - :ref:`survey_module`
     - % SURVEY


.. rubric:: Classes


.. list-table::

   * - :class:`ataddmpoleerrors`
     - ataddrandmpole adds a random multipole component to all elements of type
   * - :class:`atinsertelems`
     - ATINSERTELEMS Insert elements at given locations in a line
   * - :class:`isatlattice`
     - ISATLATTICE tests if an input argument is a valid AT lattice.
   * - :class:`setcellstruct`
     - SETCELLSTRUCT sets the field values of MATLAB cell array of structures
   * - :class:`atwritejson`
     - ATWRITEJSON Create a JSON file to store an AT lattice
   * - :class:`atGetRingProperties`
     - ATGETRINGPROPERTIES Get the ring properties
   * - :class:`atlocateparam`
     - ATLOCATEPARAM    Private function. Locate the RingParam element
   * - :class:`atsettilt`
     - ATSETTILT sets the entrance and exit rotation matrices
   * - :class:`insertelem0`
     - INSERTELEM0 - quick and dirty:
   * - :class:`atindex`
     - ATINDEX extracts the information about element families and
   * - :class:`atmaincavities`
     - ATMAINCAVITIES Get the fundamental mode cavities
   * - :class:`combinelinear45`
     - COMBINELINEAR45 combines adjacent  elements that use 4-by-5 PassMethods
   * - :class:`atwritem`
     - ATWRITEM Creates a .m file to store an AT structure
   * - :class:`atsetshift`
     - ATSETSHIFT sets the misalignment vectors
   * - :class:`findcells`
     - FINDCELLS performs a search on MATLAB cell arrays of structures
   * - :class:`atguessclass`
     - ATGUESSCLASS Tries to determine the class of an element
   * - :class:`settilt`
     - SETTILT sets the entrance and exit misalignment matrixes
   * - :class:`exitfields`
     - EXITFIELDS() Return the list of field names affecting the element exit
   * - :class:`atgetfieldvalues`
     - ATGETFIELDVALUES retrieves the field values AT cell array of elements
   * - :class:`atreduce`
     - ATREDUCE Remove useless elements from an AT structure
   * - :class:`symplectify`
     - symplectify makes a matrix more symplectic
   * - :class:`insertindrift`
     -  INSERTINDRIFT inserts one or more elements into a drift element
   * - :class:`atloadlattice`
     - ATLOADLATTICE  load a lattice from a file
   * - :class:`atfittune`
     - ATFITTUNE Fit linear tunes by scaling 2 quadrupole families
   * - :class:`atdivelem`
     - ATDIVELEM  Divide an element into pieces
   * - :class:`atsetRFCavity`
     -   ATSETRFCAVITY Set the RF Cavity with the passmethod RFCavityPass.
   * - :class:`atparamscan`
     - ATPARAMSCAN    Private function. Updates the RingParam element
   * - :class:`atSetRingProperties`
     - atSetRingProperties	Add or modify properties of the lattice
   * - :class:`atgetcells`
     - ATGETCELLS performs a search on MATLAB cell arrays of structures
   * - :class:`atsbreak`
     - ATSBREAK Insert markers at given s positions in a lattice
   * - :class:`mergedrift`
     - MERGEDRIFT removes a lattice element and merges the two adjacent drift spaces
   * - :class:`settags`
     - SETTAGS sets the 'Tag' field in AT lattice elements
   * - :class:`atsimplering`
     - 
   * - :class:`ataddmpolecomppoly`
     - ataddmpolecomppoly adds a multipole component to an existing polynomial,
   * - :class:`mvelem`
     - MVELEM Move an element
   * - :class:`setshift`
     - SETSHIFT sets the misalignment vectors T1, T2 for elements
   * - :class:`rmelem0`
     - RMELEM0 removes elements of length 0 from the accelerator lattice
   * - :class:`atelem`
     - ATELEM makes a new AT element structure from another element,
   * - :class:`mvfield`
     - MVFIELD Move fields from one structure to another
   * - :class:`findtags`
     - FINDTAGS looks for string matches in 'Tag' field of AT lattice elements
   * - :class:`buildlat`
     - BUILDLAT places elements from FAMLIST into cell array THERING
   * - :class:`atsetfieldvalues`
     - ATSETFIELDVALUES sets the field values of MATLAB cell array of structures
   * - :class:`at2str`
     - AT2STR Makes the string representation of an AT element
   * - :class:`at2py`
     - ELSTR=AT2PY(ELEM) convert AT element tp pyat
   * - :class:`combinebypassmethod`
     - COMBINEBYPASSMETHOD combines adjacent elements that have the same specified pass method
   * - :class:`atfitchrom`
     - ATFITCHROM Fit chromaticites by scaling 2 sextupole families
   * - :class:`atrotatelattice`
     - ATROTATELATTICE    Circularly shift the lattice elements
   * - :class:`attiltelem`
     - ATTILTELEM sets new rotation parameters
   * - :class:`atmakefielderrstruct`
     -  MAKERNDFIELDERRS will create a field error data structure
   * - :class:`isatelem`
     - ISATELEM tests if an input argument is a valid AT element.
   * - :class:`atwritepy`
     - ATWRITEPY Creates pyAT lattice from a Matlab lattice
   * - :class:`atsplitelem`
     - ATSPLITELEM Creates a line by inserting one or more elements into a base element
   * - :class:`splitdrift`
     - SPLITDRIFT inserts an element into a drift space
   * - :class:`atCheckRingProperties`
     - ATCHECKRINGPROPERTIES Get the ring properties if existing
   * - :class:`atshiftelem`
     - ATSHIFTELEM set new displacement parameters
   * - :class:`atfastring`
     - ATFASTRING Generate simplified AT structures
   * - :class:`entrancefields`
     - ENTRANCEFIELDS() Return the list of field names affecting the element entrance
   * - :class:`getcellstruct`
     - GETCELLSTRUCT retrieves the field values MATLAB cell array of structures
   * - :class:`atloadfielderrs`
     -  ATLOADFIELDERRS will load a field error structure into a ring

.. py:class:: atparticle


.. rubric:: Functions


.. list-table::

   * - :func:`ataddmpoleerrors`
     - ataddrandmpole adds a random multipole component to all elements of type
   * - :func:`atinsertelems`
     - ATINSERTELEMS Insert elements at given locations in a line
   * - :func:`isatlattice`
     - ISATLATTICE tests if an input argument is a valid AT lattice.
   * - :func:`setcellstruct`
     - SETCELLSTRUCT sets the field values of MATLAB cell array of structures
   * - :func:`atwritejson`
     - ATWRITEJSON Create a JSON file to store an AT lattice
   * - :func:`atGetRingProperties`
     - ATGETRINGPROPERTIES Get the ring properties
   * - :func:`atlocateparam`
     - ATLOCATEPARAM    Private function. Locate the RingParam element
   * - :func:`atsettilt`
     - ATSETTILT sets the entrance and exit rotation matrices
   * - :func:`insertelem0`
     - INSERTELEM0 - quick and dirty:
   * - :func:`atindex`
     - ATINDEX extracts the information about element families and
   * - :func:`atmaincavities`
     - ATMAINCAVITIES Get the fundamental mode cavities
   * - :func:`combinelinear45`
     - COMBINELINEAR45 combines adjacent  elements that use 4-by-5 PassMethods
   * - :func:`atwritem`
     - ATWRITEM Creates a .m file to store an AT structure
   * - :func:`atsetshift`
     - ATSETSHIFT sets the misalignment vectors
   * - :func:`findcells`
     - FINDCELLS performs a search on MATLAB cell arrays of structures
   * - :func:`atguessclass`
     - ATGUESSCLASS Tries to determine the class of an element
   * - :func:`settilt`
     - SETTILT sets the entrance and exit misalignment matrixes
   * - :func:`exitfields`
     - EXITFIELDS() Return the list of field names affecting the element exit
   * - :func:`atgetfieldvalues`
     - ATGETFIELDVALUES retrieves the field values AT cell array of elements
   * - :func:`atreduce`
     - ATREDUCE Remove useless elements from an AT structure
   * - :func:`symplectify`
     - symplectify makes a matrix more symplectic
   * - :func:`insertindrift`
     -  INSERTINDRIFT inserts one or more elements into a drift element
   * - :func:`atloadlattice`
     - ATLOADLATTICE  load a lattice from a file
   * - :func:`atfittune`
     - ATFITTUNE Fit linear tunes by scaling 2 quadrupole families
   * - :func:`atdivelem`
     - ATDIVELEM  Divide an element into pieces
   * - :func:`atsetRFCavity`
     -   ATSETRFCAVITY Set the RF Cavity with the passmethod RFCavityPass.
   * - :func:`atparamscan`
     - ATPARAMSCAN    Private function. Updates the RingParam element
   * - :func:`atSetRingProperties`
     - atSetRingProperties	Add or modify properties of the lattice
   * - :func:`atgetcells`
     - ATGETCELLS performs a search on MATLAB cell arrays of structures
   * - :func:`atsbreak`
     - ATSBREAK Insert markers at given s positions in a lattice
   * - :func:`mergedrift`
     - MERGEDRIFT removes a lattice element and merges the two adjacent drift spaces
   * - :func:`settags`
     - SETTAGS sets the 'Tag' field in AT lattice elements
   * - :func:`atsimplering`
     - 
   * - :func:`ataddmpolecomppoly`
     - ataddmpolecomppoly adds a multipole component to an existing polynomial,
   * - :func:`mvelem`
     - MVELEM Move an element
   * - :func:`setshift`
     - SETSHIFT sets the misalignment vectors T1, T2 for elements
   * - :func:`rmelem0`
     - RMELEM0 removes elements of length 0 from the accelerator lattice
   * - :func:`atelem`
     - ATELEM makes a new AT element structure from another element,
   * - :func:`mvfield`
     - MVFIELD Move fields from one structure to another
   * - :func:`findtags`
     - FINDTAGS looks for string matches in 'Tag' field of AT lattice elements
   * - :func:`buildlat`
     - BUILDLAT places elements from FAMLIST into cell array THERING
   * - :func:`atsetfieldvalues`
     - ATSETFIELDVALUES sets the field values of MATLAB cell array of structures
   * - :func:`at2str`
     - AT2STR Makes the string representation of an AT element
   * - :func:`at2py`
     - ELSTR=AT2PY(ELEM) convert AT element tp pyat
   * - :func:`combinebypassmethod`
     - COMBINEBYPASSMETHOD combines adjacent elements that have the same specified pass method
   * - :func:`atfitchrom`
     - ATFITCHROM Fit chromaticites by scaling 2 sextupole families
   * - :func:`atrotatelattice`
     - ATROTATELATTICE    Circularly shift the lattice elements
   * - :func:`attiltelem`
     - ATTILTELEM sets new rotation parameters
   * - :func:`atmakefielderrstruct`
     -  MAKERNDFIELDERRS will create a field error data structure
   * - :func:`isatelem`
     - ISATELEM tests if an input argument is a valid AT element.
   * - :func:`atwritepy`
     - ATWRITEPY Creates pyAT lattice from a Matlab lattice
   * - :func:`atsplitelem`
     - ATSPLITELEM Creates a line by inserting one or more elements into a base element
   * - :func:`splitdrift`
     - SPLITDRIFT inserts an element into a drift space
   * - :func:`atCheckRingProperties`
     - ATCHECKRINGPROPERTIES Get the ring properties if existing
   * - :func:`atshiftelem`
     - ATSHIFTELEM set new displacement parameters
   * - :func:`atfastring`
     - ATFASTRING Generate simplified AT structures
   * - :func:`entrancefields`
     - ENTRANCEFIELDS() Return the list of field names affecting the element entrance
   * - :func:`getcellstruct`
     - GETCELLSTRUCT retrieves the field values MATLAB cell array of structures
   * - :func:`atloadfielderrs`
     -  ATLOADFIELDERRS will load a field error structure into a ring

.. py:function:: ataddmpoleerrors

   | ataddrandmpole adds a random multipole component to all elements of type
   | 'type' where type can be 'dipole', 'quadrupole', or 'sextupole'
   | 
   | [newring] = ATRANDMPOLE(ring,type,newindex,strength,radius)
   | 
   | ring = input ring
   | type = 'dipole', 'quadrupole' or 'sextupole'
   | newindex: index of Multipole to add
   | strength: strength of the multipole component at the given radius
   | radius: reference radius for defining the absolute strength
   | if randflag is set to 1, then the errors will be random, Gaussian
   | distributed
   |  The formula for the added errors is
   |  B^(N)_(n) = radius^(n-N)*b_n/b_N
   |  It represents the relative field error to the design value at the ref.
   |  radius
   | For example, to add a random octupole error of 1e-4 at 25 mm, relative to all
   | quadrupoles:
   |  newring =ataddrandmpole(ring,'quadrupole',4,.1e-4,.025,1);
   | 
   | See also: ataddmpolecomppoly attiltelem atshiftelem

.. py:function:: atinsertelems

   | ATINSERTELEMS Insert elements at given locations in a line
   | 
   | NEWLINE=ATINSERTELEMS(LINE,REFPTS,FRAC1,ELEM1[,FRAC2,ELEM2...])
   |    a new line is created by inserting elements at each location specified
   |    by REFPTS.
   | 
   | LINE:      Cell array of structures
   | REFPTS:    Insertion points (index list or logical mask)
   | FRAC:      Location of the inserted element ELEM within LINE{REFPTS(i)}
   |            0<=FRAC<=1
   |  if FRAC = 0, ELEM is inserted before LINE{REFPTS(i)} (no splitting)
   |  if FRAC = 1, ELEM is inserted after LINE{REFPTS(i)} (no splitting)
   |  if FRAC = NaN, LINE{REFPTS(i)} is replaced by ELEM (no check for identical length)
   |  if ELEM = [], nothing is inserted, only the splitting takes place
   | 
   |  FRAC and ELEM must be scalars or array of the same size as REFPTS
   | 
   |  See also ATSPLITELEM ATDIVELEM

.. py:function:: isatlattice

   | ISATLATTICE tests if an input argument is a valid AT lattice.
   | 
   |   A valid AT lattice is a MATLAB cell array of valid AT elements
   | 
   |   [TEST, BADLIST, ERRORSTR] = ISATLATTICE(ELEM, ['option1',...])
   | 
   |   Allowed otions:
   |    'display' - print problems to command window;
   |    'stop'    - return after the first problem is found
   | 
   |   TEST     - test result,  1 = valid AT element
   |   ERRORSTR - multi-line error message
   | 
   |   See also ISATELEM, ATELEM, ATLATTICE

.. py:function:: setcellstruct

   | SETCELLSTRUCT sets the field values of MATLAB cell array of structures
   | 
   |    Note that the calling syntax must be in the form of assignment:
   |    CELLARRAY = SETCELLSTRUCT(CELLARRAY,...)
   |    MATLAB does not modify variables that only appear on the right
   |    hand side as arguments.
   | 
   |  Numeric data
   |  ---------------------------------------------------------
   |  CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE,M,N)
   | 	 Sets (M,N) element equal to VALUE when the field data is
   |    a matrix. The assigned VALUE may be
   |    1. Scalar numeric value to be written to all CELLARRAY elements
   |    2. Numeric array of the same length as INDEX array
   | 
   |  CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE,M) can be used
   |    for one dimensional vectors
   | 
   |  CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE) is the same as
   |    SETCELLSTRUCT(CELLARRAY,'field',index,1,1) if the field data
   |    is a scalar
   | 
   |  Character array
   |  --------------------------------------------------------------------
   |  CELLARRAY SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE) is a MATLAB
   |    cell array of strings when specified fields contain strings.
   |    The assignment VALUE may be
   |    1. Character string,
   |    2. Character array with the number of rows matching the number of
   |        elements in INDEX array,
   |    3. Cell array of strings, with either one element or with the same
   |        length as index.
   | 
   |  See also ATSETFIELDVALUES GETCELLSTRUCT FINDCELLS

.. py:function:: atwritejson

   | ATWRITEJSON Create a JSON file to store an AT lattice
   | 
   | JS=ATWRITEJSON(RING)
   |    Return the JSON representation of RING as a character array
   | 
   | ATWRITEJSON(RING, FILENAME)
   |    Write the JSON representation of RING to the file FILENAME
   | 
   | ATWRITEJSON(RING, ..., 'compact', true)
   |    If compact is true, write a compact JSON file (no linefeeds)
   | 
   | see also atloadlattice

.. py:function:: atGetRingProperties

   | ATGETRINGPROPERTIES Get the ring properties
   | 
   |  [V1,V2,...]=ATGETRINGPROPERTIES(RING,'param1','param2',...)
   |    Extract lattice properties. Extract from the RingParam element of the
   |    lattice if present, or from the lattice elements. The number of outputs
   |    corresponds to the number of property names in input.
   | 
   |  RING:             Ring structure
   | 
   |  Standard properties (case independent names):
   |    'FamName'               Name of the lattice
   |    'name'                   "   "   "     "
   |    'Energy'                Ring energy [eV]
   |    'Periodicity'           Number of periods to build a full ring
   |    'Particle'              particle (Particle object)
   |    'cavpts'                Location of the main cavities
   |    'beta'                  Relativistic beta of the particles
   |    'gamma'                 Relativistic gamma of the particles
   |    'rf_frequency'          RF frequency (main cavities) [Hz]
   |    'rf_timelag'            RF timelag (main cavities) [m]
   |    'BRho'                  Particle rigidity [T.m]
   |    'mcf'                   Momentum compaction factor "alpha"
   |    'slip_factor'           Slip factor "eta"
   |    'is_6d'                 Longitudinal motion (cavities, radiating elements, ...)
   |    'radiation'                  "          "
   |    'has_cavity'            Presence of an active RF cavity (implies 'radiation')
   | 
   |  Properties for the full ring ('periods' x cells):
   |    'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
   |    'harmonic_number'          "       "
   |    'Circumference'         Ring circumference [m] (cell_length * Periodicity)
   |    'rf_voltage'            RF voltage [V] (cell_rf_voltage * Periodicity)
   |    'revolution_frequency'  Revolution frequency [Hz] (cell_revolution_frequency / Periodicity)
   | 
   |  Properties for one cell:
   |    'cell_harmnumber'       Harmonic number (for 1 cell)
   |    'cell_length'           Cell length [m]
   |    'cell_rf_voltage'       RF voltage per cell [V] (main cavities)
   |    'cell_revolution_frequency' Revolution frequency [Hz] (for 1 cell)
   | 
   |  Custom properties may be added with atSetRingProperties. Custom property
   |  names are case dependent.
   | 
   |  Example:
   |  >> [energy, beta] = atGetRingProperties(ring,'Energy','beta');
   | 
   |  PROPERTIES=ATGETRINGPROPERTIES(RING)
   | 
   |  RING:             Ring structure
   |  PROPERTIES:       Structure with fields:
   |    'FamName'               Name of the lattice
   |    'Energy'                Ring energy [eV]
   |    'Periodicity'           Number of periods to build a full ring
   |    'Particle'              particle (Particle object)
   |    'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
   |    'cell_harmnumber'       Harmonic number (for 1 cell)
   |    'cavpts'                Location of the main cavities
   | 
   |  PROPERTIES=ATGETRINGPROPERTIES(RING,'all')
   |    The PROPERTIES structure contains all the properties described above,
   |    both standard and custom.
   | 
   |  For fast access, some ring properties are stored in a RingParam element
   |  ideally located in the 1st position of the lattice. Without such element,
   |  the properties are deduced from the lattice contents. This is much slower
   |  and ATGETRINGPROPERTIES displays a warning indicating how to add the
   |  RingParam element:
   | 
   | >>ring=atSetRingProperties(ring)
   | 
   |   See also ATSETRINGPROPERTIES

.. py:function:: atlocateparam

   | ATLOCATEPARAM    Private function. Locate the RingParam element
   | 
   |  IDX=ATLOCATEPARAM(ring)
   | 
   |  ring:         Ring structure
   |  IDX:          Index of the 1st RingParam element in the ring
   | 
   |   See also ATGETRINGPROPERTIES, ATSETRINGPROPERTIES

.. py:function:: atsettilt

   | ATSETTILT sets the entrance and exit rotation matrices
   |  of an element or a group of elements in THERING
   | 
   |  RING=ATSETTILT(RING,ELEMINDEX, PSI)
   |  ELEMINDEX contains indexes of elements to be rotated
   |  PSI - angle(s) of rotation in RADIANS
   |    POSITIVE PSI corresponds to a CORKSCREW (right)
   |    rotation of the ELEMENT looking in the direction of the beam.
   |    (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
   |    The misalgnment matrixes are stored in fields R1 and R2
   |    R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
   |    R2 = R1'
   | 
   |  RING=ATSETTILT(RING,ELEMINDEX,PSI,'RelativeTilt')
   |  the rotation is added to the previous misalignment
   | 
   |  ATSETTILT(ELEMINDEX, PSI) Uses the global variable THERING
   | 
   |  See also ATSETSHIFT

.. py:function:: insertelem0

   | INSERTELEM0 - quick and dirty:
   |   inserts element(s) of zero length into a drift space
   | 
   |  NEWLATTICE = INSERTELEM0(LATTICE, DRIFTINDEX, SPLITLENGTH, ELEMDATA)

.. py:function:: atindex

   | ATINDEX extracts the information about element families and
   |  indexes from AT lattice
   | 
   | ATI=ATINDEX(LATTICE)
   |   returns a srtucture with fields named after element family
   |   containing an array of their AT indexes;
   | 
   |    ATI.QF = [...]
   |    ATI.QD = [...];
   |    ...
   | 
   | ATI=ATINDEX() Uses the global variable THERING
   | 
   |  See also ATGETCELLS

.. py:function:: atmaincavities

   | ATMAINCAVITIES Get the fundamental mode cavities
   | 
   | IDCAV=ATMAINCAVITIES(RING)
   |    Return a refpts pointing to fundamental mode cavities
   |    (cavities with the lowest frequency)

.. py:function:: combinelinear45

   | COMBINELINEAR45 combines adjacent  elements that use 4-by-5 PassMethods
   |  [NEWLATTICE, SHIFTEDKEEPINDEX, SHIFTEDREF] = COMBINELINEAR45(LATTICE,KEEPINDEX,REFINDEX)

.. py:function:: atwritem

   | ATWRITEM Creates a .m file to store an AT structure
   | 
   |   INPUTS
   |   1. ring    - Lattice structure (.mat file)
   |   2. filname - output file where to write the lattice as a line by line file using
   |   element constructors
   | 
   |   EXAMPLES
   |   1. atwritem(ring) % Prints the result in the command window
   |   2. atwritem(ring,filename) % Prints the result in a file
   | 
   |   See also at2str

.. py:function:: atsetshift

   | ATSETSHIFT sets the misalignment vectors
   | 
   |  RING=ATSETSHIFT(RING, ELEMINDEX, DX, DY) sets the entrance and exit misalignment vectors
   |   of one element or a group of elements in the globally defined lattice THERING.
   | 
   |  ELEMINDEX contains indexes of elements to be misaligned
   |  DX, DY are displacements of the ELEMENT
   |   so the coordinate transformation on the particle at entrance is
   | 	X  ->  X-DX
   |    Y  ->  Y-DY
   | 
   |  RING=ATSETSHIFT(RING, ELEMINDEX, DX, DY, 'RelativeShift')
   |  adds the shifts DX and DY to the the previous misalignments
   |  of the elements
   | 
   |  ATSETSHIFT(ELEMINDEX, DX, DY) Uses the global variable THERING
   | 
   |  See also: ATSETTILT

.. py:function:: findcells

   | FINDCELLS performs a search on MATLAB cell arrays of structures
   | 
   |  INDEX = FINDCELLS(CELLARRAY, 'field')
   |    returns indexes of elements that have a field named 'field'
   | 
   |  INDEX = FINDCELLS(CELLARRAY, 'field', VALUE)
   |    returns indexes of elements whose field 'field'
   |    is equal to VALUE1, VALUE2, ... or VALUEN. Where VALUE can either be
   |    character strings or a number. If its a character string REGULAR
   |    expressions can be used.
   | 
   |  Example:
   |    findcells(THERING,'Length',0, 0.2);  % will match elements of
   |                                           lengths 0 and 0.2
   |    findcells(THERING,'FamName','SFA','SDA');
   | 
   |  See also ATGETCELLS, GETCELLSTRUCT, SETCELLSTRUCT, REGEXPI

.. py:function:: atguessclass

   | ATGUESSCLASS Tries to determine the class of an element
   | ATCLASS=ATGUESSCLASS(ATELEM) Tries to determine the class of an element
   | 
   |   INPUTS
   |   1. elem - AT element
   | 
   | 
   |   NOTES
   |   1. atclass=atguessclass(atelem,'useclass')
   |   By default, ATGUESSCLASS will default "useless" elements (PolynopmB==0)
   |   to 'Drift' or 'Marker', depending on 'Length'. When specifying
   |   'UseClass', it it will preserve the 'Class' field for those elements.
   | 
   |   See also atwritem

.. py:function:: settilt

   | SETTILT sets the entrance and exit misalignment matrixes
   |  of an element or a group of elements in THERING
   |  Previously stored values are overwritten.
   | 
   |  SETTILT(ELEMINDEX, PSI)
   |  ELEMINDEX contains indexes of elements to be rotated
   |  PSI - angle(s) of rotation in RADIANS
   |    POSITIVE PSI corresponds to a CORKSCREW (right)
   |    rotation of the ELEMENT.
   |    (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
   |    The misalgnment matrixes are stored in fields R1 and R2
   |    R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
   |    R2 = R1'
   | 
   |   NOTES
   |   1. This function is deprecated. Use atsettilt instead
   | 
   |  See also SETSHIFT MKSROLLMAT

.. py:function:: exitfields

   | EXITFIELDS() Return the list of field names affecting the element exit
   | 
   |  Optional arguments:
   |  'KeepAxis', if present, rotations translations are excluded from list
   | 
   | see also: entrancefields atdivelem

.. py:function:: atgetfieldvalues

   | ATGETFIELDVALUES retrieves the field values AT cell array of elements
   | 
   |  VALUES = ATGETFIELDVALUES(RING,'field') extracts the values of
   |  the field 'field' in all the elements of RING
   | 
   |  VALUES = ATGETFIELDVALUES(RING,INDEX,'field') extracts the values of
   |  the field 'field' in the elements of RING selected by INDEX
   | 
   |  if RING{I}.FIELD is a numeric scalar
   |     VALUES is a length(INDEX) x 1 array
   |  otherwise
   |     VALUES is a length(INDEX) x 1 cell array
   | 
   | VALUES = ATGETFIELDVALUES(...,'Default',default_value) Uses default_values
   |    if the required field is not existing
   | 
   |  More generally ATGETFIELDVALUES(RING,INDEX,subs1,subs2,...) will call
   |   GETFIELD(RING{I},subs1,subs2,...) for I in INDEX
   | 
   |  Examples:
   | 
   |  V=ATGETFIELDVALUES(RING,1:10,'PolynomB') is a 10x1  cell array
   |  such that V{I}=RING{I}.PolynomB for I=1:10
   | 
   |  V=ATGETFIELDVALUES(RING(1:10),'PolynomB',{1,2}) is a 10x1 array
   |  such that V(I)=RING{I},PolynomB(1,2)
   | 
   | 
   |  See also ATSETFIELDVALUES ATGETCELLS GETCELLSTRUCT FINDCELLS

.. py:function:: atreduce

   | ATREDUCE Remove useless elements from an AT structure
   | NEWRING=ATREDUCE(RING)
   | 
   |  Remove elements with PassMethod='IdentityPass' and merges adjacent
   |  similar elements
   | 
   | NEWRING=ATREDUCE(RING,REFPTS)
   | 	When merging similar elements, keep REFPTS intact.
   | 
   | [NEWRING,KEPT]=ATREDUCE(...)
   | 	Returns the index of kept elements so that NEWRING=OLDRING(KEPT)
   | 
   | [NEWRING,KEPT,NEWREFPTS]=ATREDUCE(RING,REFPTS)
   | 	Returns in addition the updated list of reference points.
   | 

.. py:function:: symplectify

   | symplectify makes a matrix more symplectic
   | follow Healy algorithm as described by McKay
   | BNL-75461-2006-CP

.. py:function:: insertindrift

   |  INSERTINDRIFT inserts one or more elements into a drift element
   |   and returns a sequence (cell array) of elements  ready to to be used
   |   in AT lattice
   | 
   |  ELEMSEQ = INSERTELEM(DRIFT0, ELEM1, POS1, ... ELEMN, POSN)
   | 
   |  EXAMPLE: FODO cell
   | 
   |  --- 1. Declare elements
   | 
   |  D  = atelem('drift','Length',4.5);
   |  QF = atelem('quad','Length', 1, 'K',  1.234);
   |  QD = atelem('quad','Length', 1, 'K', -2.345);
   | 
   |  --- 2. Insert quadrupoles in the drift;
   | 
   |  FODOCELL = insertindrift(D, QF, 0.5, QD, 2, QF, 3.5);
   | 
   |  See also: SPLITELEM

.. py:function:: atloadlattice

   | ATLOADLATTICE  load a lattice from a file
   | 
   | LATTICE=ATLOADLATTICE(FILEPATH)
   |    Create an AT lattice (cell array of AT elements) from FILEPATH
   | 
   | LATTICE=ATLOADLATTICE(FILEPATH,'matkey',VARIABLE_NAME,...)
   |    Give the name of a variable containing the lattice in .mat files
   |    containing several variables. By default, if the .mat file contains a
   |    single variable, it will be loaded. Otherwise ATLOADLATTICE will look
   |    in order for 'ring', 'lattice' and 'RING'.
   | 
   | LATTICE=ATLOADLATTICE(FILEPATH,KEY,VALUE...)
   |    (Key,value) pairs are added to the RingProperties of the lattice or
   |    overload existing ones.
   |    Standard keys include: FamName, Energy, Periodicity, HarmNumber, Particle
   |    Custom properties are allowed
   | 
   | Allowed file types:
   |    .mat    Binary Matlab file. If the file contains several variables,
   |            the variable name must be specified using the 'matkey' keyword.
   | 
   |    .m      Matlab function. The function must output a valid AT structure.
   |    .json   JSON file
   | 
   | see also atwritem, atwritejson

.. py:function:: atfittune

   | ATFITTUNE Fit linear tunes by scaling 2 quadrupole families
   | 
   |  NEWRING = ATFITTUNE(RING,NEWTUNES,QUADFAMILY1,QUADFAMILY2)
   |  NEWRING = ATFITTUNE(RING,DPP,NEWTUNES,QUADFAMILY1,QUADFAMILY2)
   | 
   | RING:          Cell array
   | DPP:           Optional momentum deviation (default 0)
   | NEWTUNES:      Desired tune values
   | QUADFAMILY1:   1st quadrupole family
   | QUADFAMILY2:   2nd quadrupole family
   | 
   | QUADFAMILY may be:
   |    string: Family name
   |    logical array: mask of selected elements in RING
   |    Numeric array: list of selected elements in RING
   |    Cell array: All elements selected by each cell
   | 
   | NEWRING = ATFITTUNE(RING,...,'UseIntegerPart') With this flag, the
   |    function fits the tunes to the total values of NEWTUNES, including
   |    the integer part.
   |    With this option the function is substantially slower!
   | 
   | NEWRING = ATFITTUNE(RING,...,'KStep',kstep)
   |    kstep is the quadrupole strength applied to build the jacobian [m^-2].
   |    Default: 1.0e-6
   | 
   | NEWRING = ATFITTUNE(RING,...,'dp',DP)
   |    Specify off-momentum if radiation is OFF (default 0)
   | 
   | NEWRING = ATFITTUNE(RING,...,'dct',DCT)
   |    Specify path lengthening if radiation is OFF (default 0)
   | 
   |  See also ATFITCHROM

.. py:function:: atdivelem

   | ATDIVELEM  Divide an element into pieces
   | 
   | LINE=ATDIVELEM(ELEM,FRAC)
   | 
   | ELEM:  Element to be divided
   | FRAC:  Length of each piece, as a fraction of the total length
   | 
   | LINE:  Cell array
   | 
   |  The sum of frac elements does not need to be 1, however for bending
   |  magnets, the length will be divided according to FRAC, but the total
   |  bending angle will be divided according to FRAC/SUM(FRAC) so that the
   |  total bending angle is kept.
   | 
   |  Example:
   | 
   | >> qf=atquadrupole('QF',0.1,0.5);
   | >> line=atdivelem(qf,[0.5;0.5]); % Split a quadrupole in two halves
   | 
   |  Optional arguments:
   |  'KeepAxis', if present, rotations translations are kept at all slices
   | 
   |  See also ATINSERTELEMS ATSLICE ATSPLITELEM ENTRANCEFIELDS EXITFIELDS

.. py:function:: atsetRFCavity

   |   ATSETRFCAVITY Set the RF Cavity with the passmethod RFCavityPass.
   |   RFCavityPass allows to change the energy of the beam changing the
   |   frequency of the cavity.
   | 
   |    newring = atSetRFCavity(ring, rfv, radflag, HarmNumber, DeltaFreq)
   |    sets the synchronous phase of the cavity, the voltage, the harmonic
   |    number and the PassMethod.
   | 
   |   INPUTS
   |   1. ring       - Ring structure
   |   2. rfv        - RF-voltage in volts
   |   3. radflag    - 0/1: activat/desactivate radiation (atradon/atradoff)
   |   4. HarmNumber - Harmonic number
   |   5. DeltaFreq  - Frequency shift in Hz
   | 
   |   OUTPUTS
   |   1. newring - update ring structure with nw RF parameters
   | 
   |   NOTES
   |   1. All the N cavities will have a voltage rfv/N
   |   2. radflag says whether or not we want radiation on, which affects
   |      synchronous phase. If radflag is 0, the function calls atradoff,
   |      if it is 1, it calls atradon.
   |   3. Cavities in the ring must have the Class RFCavity.
   |   4. Normally DeltaFreq should be 0, it's different from 0 when you want to
   |      simulate a change of energy changing the RF frequency. DeltaFreq is in
   |       Hz.
   |   5. Does not work well for misaligned cavity
   | 
   |   EXAMPLES
   | 
   |    1. normal use:
   |    newring = atsetRFCavity( ring, 6e6, 1, 992, 0 )
   | 
   |    2. for off-energy simulations:
   |    newring = atsetRFCavity( ring, 6e6, 1, 992, 100 )
   | 
   |   See also atsetcavity, atradon, atradoff, atgetU0

.. py:function:: atparamscan

   | ATPARAMSCAN    Private function. Updates the RingParam element
   | 
   | NEWPARMS=ATPARAMSCAN(RING,PARMS,VARARGIN)
   | 
   |   See also ATGETRINGPROPERTIES, ATSETRINGPROPERTIES

.. py:function:: atSetRingProperties

   | atSetRingProperties	Add or modify properties of the lattice
   | 
   |  newring=atSetRingProperties(ring [,key,value]...)
   |    Add or modify the attributes of the RingParam element of the lattice,
   |    Insert a new RingParam element if necessary
   | 
   |  Available properties (property names are case-independent):
   |    'FamName'               Name of the lattice
   |    'name'                   "   "   "     "
   |    'Energy'                Ring energy [eV]
   |    'Periodicity'           Number of periods to build a full ring
   |    'Particle'              particle (perticle name or Particle object)
   |    'cavpts'                Location of the main cavities
   |    'rf_frequency'          RF frequency (main cavities) [Hz]. Use 'nominal'
   |                            to set the nominal frequency
   |    'rf_timelag'            RF timelag (main cavities) [m]
   | 
   |  Properties for the full ring ('periods' x cells):
   |    'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
   |    'harmonic_number'          "       "
   |    'rf_voltage'            RF voltage [V] (cell_rf_voltage * Periodicity)
   | 
   |  Properties for one cell:
   |    'cell_harmnumber'       Harmonic number (for 1 cell)
   |    'cell_rf_voltage'       RF voltage per cell [V] (main cavities)
   | 
   |  Additional custom fields may be added. They can be retrieved by
   |  atGetRingProperties and are saved in files. Custom field names are
   |  case-dependent.
   | 
   |  For fast access, some ring properties are stored in a RingParam element
   |  ideally located in the 1st position of the lattice. If there is no such
   |  element, atSetRingProperties will add it.
   | 
   |  See also atGetRingProperties

.. py:function:: atgetcells

   | ATGETCELLS performs a search on MATLAB cell arrays of structures
   | 
   |  OK = ATGETCELLS(RING, 'field')
   |    returns indexes of elements that have a field named 'field'
   | 
   |  OK = ATGETCELLS(RING, 'field', VALUE1...)
   |    returns indexes of elements whose field 'field'
   |    is equal to VALUE1, VALUE2, ... or VALUEN. Where VALUE can either be
   |    character strings or a number. If its a character string REGULAR
   |    expressions can be used.
   | 
   |  OK = ATGETCELLS(RING, 'field', @TESTFUNCTION, ARGS...)
   |    Uses the user-defined TESTFUNCTION to select array elements
   |    TESTFUNCTION must be of the form:
   |        OK=TESTFUNTION(ATELEM,FIELDVALUE,ARGS...)
   | 
   |  OK is a logical array with the same size as RING, refering to matching
   |  elements in RING
   | 
   |  See also ATGETFIELDVALUES, ATSETFIELDVALUES, FINDCELLS, REGEXPI

.. py:function:: atsbreak

   | ATSBREAK Insert markers at given s positions in a lattice
   | 
   | [NEWRING,REFPTS]=ATSBREAK(RING,SPOS) build a new lattice with breakpoints
   | 
   | RING:  AT structure
   | SPOS:  Poition of the desired breakpoints
   | 
   | NEWRING:   Modified AT structure
   | REFPTS:    Index of breakpoints
   | 
   |  See also ATINSERTELEMS ATSPLITELEM ATDIVELEM

.. py:function:: mergedrift

   | MERGEDRIFT removes a lattice element and merges the two adjacent drift spaces
   | 
   |  MERGEDRIFT (SPLITPOS) removes an element located at SPLITPOS from the global lattice THERING
   |  surrounded by two DRIFT spaces. The resulting drift has Length L0 = L1 + LSPLIT + L2;
   |  Number of elements in THERING is thus reduced by 2
   | 
   |  See also: SPLITDRIFT

.. py:function:: settags

   | SETTAGS sets the 'Tag' field in AT lattice elements
   |  LATTICE = SETTAGS(LATTICE, INDEX, TAG)
   |    INDEX can be integer AT index or a string famly name
   |    TAG is a string tag or a cell array of strings
   |  LATTICE = SETTAGS(LATTICE, INDEX, TAG, 'append')
   |    appends to existing tags

.. py:function:: atsimplering


.. py:function:: ataddmpolecomppoly

   | ataddmpolecomppoly adds a multipole component to an existing polynomial,
   | scaling it by the reference multipole value and a radius to set the length
   | scale
   | [PolynomOut] = ATADDMPOLECOMPPOLY(Polynom,index,strength,radius)
   | 
   | Polynom = given field polynomial that we want to add a component to
   | refindex = index of magnet: 1 for dipole, 2 for quadrupole
   | newindex: index of Multipole to add
   | strength: strength of the multipole component at the given radius
   | radius: reference radius for defining the absolute strength
   |  B^(N)_(n) = radius^(n-N)*b_n/b_N
   | 
   | See also: attiltelem, atshiftelem, atmodelem

.. py:function:: mvelem

   | MVELEM Move an element
   | 
   | MVELEM(ELEMPOS, DIST) Move an element located at ELEMPOS in THERING
   |  surrounded by two DRIFT spaces
   | 
   |   0   < DIST  < LD move downstream
   |  -LU  < DIST  < 0  move upstream
   |   where LU and LD - lenths of
   |   upstream and downstrem drift drifts BEFORE!!! the move
   | 
   |  Number of elements in THERING and total length remain the same
   | 
   |  See also: SPLITDRIFT, MERGEDRIFT

.. py:function:: setshift

   | SETSHIFT sets the misalignment vectors T1, T2 for elements
   | 
   |  SETSHIFT(ELEMINDEX, DX, DY) sets the entrance and exit misalignment vectors
   |   of one element or a group of elements in the globally defined lattice THERING.
   | 
   |   DX, DY are displacements of the ELEMENT
   |   so the coordinate transformation on the particle at entrance is
   | 	X  ->  X-DX
   |    Y  ->  Y-DY
   |   The elements to be modified are given by ELEMINDEX
   | 	Previous stored values are overwritten.
   | 
   |  See also SETSHIFT

.. py:function:: rmelem0

   | RMELEM0 removes elements of length 0 from the accelerator lattice
   |  NEWLATTICE = RMELEM0(LATTICE,ELEMINDEX)
   |  [NEWLATTICE, SHIFTEDINDEX] = RMELEM0(LATTICE,ELEMINDEX)
   | 
   |  The number of elements in the modified lattice is
   |  reduced by length(ELEMINDEX)
   | 
   |  SHIFTEDINDEX points to elements in the NEWLATTICE that
   |  immediateley followed the removed elements in the original LATTICE
   | 
   |  See also: SPLITDRIFT MERGEDRIFT

.. py:function:: atelem

   | ATELEM makes a new AT element structure from another element,
   |  standard AT type, or a template on file. Existing fields are overridden
   |  with new values. New fields are added and initialized.
   | 
   |  NEWELEM = ATELEM(ELEM,'Field1','Value1','Field2', 'Value2', ...)
   |   ELEM can be 1) another element structure
   |               2) One of the standard AT types
   |                  'drift'
   |                  'quadrupole'
   |                  'sextupole'
   |                  'bend','rbend','sbend'
   |                  'marker'
   |                  'corrector'
   |                   ...

.. py:function:: mvfield

   | MVFIELD Move fields from one structure to another
   | 
   | [NEWDST,NEWSRC]=MVFIELD(DST,SRC,FIELDNAMES)
   |    Moves the selected fields from SRC to DST
   | 
   | DST:           Destination structure
   | SRC:           Source structure
   | FIELDNAMES:    Field names to be moved (Default: all fields from SRC)
   | 
   | NEWDST:        DST structure with fields added
   | NEWSRC:        SRC structure with fields removed

.. py:function:: findtags

   | FINDTAGS looks for string matches in 'Tag' field of AT lattice elements
   | 
   |  INDEX = FINDTAGS(CELLARRAY, MATCHSTR)
   |    returns indexes of elements that have a field 'Tag'
   |    whose value is a string exactly matching MATCHSTR
   |    or a cell array of strings with one element matching MATCHSTR
   | 
   |  See also FINDCELLS, SETTAGS,

.. py:function:: buildlat

   | BUILDLAT places elements from FAMLIST into cell array THERING
   |  in the order given by ineger arry ELIST
   |  to be use in Accelerator Toolbox lattice definition files

.. py:function:: atsetfieldvalues

   | ATSETFIELDVALUES sets the field values of MATLAB cell array of structures
   | 
   |    Note that the calling syntax must be in the form of assignment:
   |    RING = ATSETFIELDVALUES(RING,...)
   |    MATLAB does not modify variables that only appear on the right
   |    hand side as arguments.
   | 
   | NEWRING=ATSETFIELDVALUES(RING,'field',VALUES)
   |    In this mode, the function will set values on all the elements of RING
   | 
   | NEWRING=ATSETFIELDVALUES(RING,INDEX,'field',VALUES)
   |    In this mode, the function will set values on the elements of RING
   |    specified by INDEX, given as a list of indices or as a logical mask
   | 
   | NEWRING=ATSETFIELDVALUES(RING,'field',VALUESTRUCT)
   |    In this mode, the function will set values on the elements of RING
   |    whose family names are given by the field names of VALUESTRUCT
   | 
   | NEWRING=ATSETFIELDVALUES(RING,RINGINDEX,....,VALUESTRUCT)
   |    As in the previous mode, the function will set values on the elements
   |    of RING whose family names are given by the field names of VALUESTRUCT.
   |    But RINGINDEX=atindex(RING) is provided to avoid multiple computations.
   | 
   |  Field selection
   |  ---------------------------------------------------------
   |  NEWRING = ATSETFIELD(RING,'field',VALUES)
   |    For each I=1:length(RING), set RING{I}.FIELD=value
   | 
   |  NEWRING = ATSETFIELD(RING,'field',{M,N},VALUES)
   |    For each I=1:length(RING), set RING{I}.FIELD(M,N)=value
   | 
   |  More generally,
   |  NEWRING = ATSETFIELD(RING,subs1,subs2,...,VALUES)
   |    For each I=1:length(RING), SETFIELD(RING{I},subs1,subs2,...,value)
   | 
   |  The first dimension of VALUES must be either length(INDEX) or 1 (the value
   |  will be repeated for each element). For a vector to be repeated, enclose
   |  it in a cell array.
   | 
   |  Value format
   |  ---------------------------------------------------------
   |  Cell array VALUES
   |  -----------------
   |  Mx1 or 1xM cell array : one cell per element
   |  1x1 cell array : cell 1 is affected to all selected elements
   | 
   |  Character array VALUES
   |  ---------------------
   |  1xN char array (string) : the string as affected to all selected elements
   |  MxN char array : one row per element
   | 
   |  Numeric array VALUES
   |  --------------------
   |  1x1 (scalar) : the value is affected to all selected elements
   |  Mx1 (column) : one value per element
   |  MxN (matrix) : one row affected per element. If M==1, the single row
   |                 is affected to all selected elements
   |  To assign column vectors to parameters, use a cell array VALUES, with one
   |  column per cell
   | 
   |  See also ATGETFIELDVALUES ATGETCELLS SETCELLSTRUCT FINDCELLS

.. py:function:: at2str

   | AT2STR Makes the string representation of an AT element
   | 
   | AT2STR Creates a string such that EVAL(AT2STR(ELEM)) recreates an
   | identical element.
   | 
   |   INPUTS
   |   1. elem - Elem to write
   | 
   |   OUTPUTS
   |   1. elsstr - String given the AT constructor of the element
   | 
   |   See also atwritem atwritepy

.. py:function:: at2py

   | ELSTR=AT2PY(ELEM) convert AT element tp pyat
   | 
   | AT2PY Creates a pyat element

.. py:function:: combinebypassmethod

   | COMBINEBYPASSMETHOD combines adjacent elements that have the same specified pass method
   |  [NEWLATTICE, SHIFTEDKEEPINDEX, SHIFTEDREF] = COMBINEBYPASSMETHOD(LATTICE,METHOD,KEEPINDEX,REFINDEX)

.. py:function:: atfitchrom

   | ATFITCHROM Fit chromaticites by scaling 2 sextupole families
   | 
   |  NEWRING = ATFITCHROM(RING,NEWCHROM,SEXTFAMILY1,SEXTFAMILY2)
   |  NEWRING = ATFITCHROM(RING,DPP,NEWCHROM,SEXTFAMILY1,SEXTFAMILY2)
   | 
   | RING:          Cell array
   | DPP:           Optional momentum deviation (default 0)
   | NEWCHROM:      Desired non-normalized chromaticities
   | SEXTFAMILY1:   1st sextupole family
   | SEXTFAMILY2:   2nd sextupole family
   | 
   | SEXTFAMILY may be:
   |    string: Family name
   |    logical array: mask of selected elements in RING
   |    Numeric array: list of selected elements in RING
   |    Cell array: All elements selected by each cell
   | 
   | NEWRING = ATFITCHROM(RING,...,'DPStep',dpstep)
   |    dpstep is the momentum step applied to compute the chromaticity.
   |    The default is the DPStep global option, which defaults to 3.0e-6
   | 
   | NEWRING = ATFITCHROM(RING,...,'HStep',hstep)
   |    hstep is the sextupole strength applied to build the jacobian [m^-3].
   |    Default: 1.0e-5
   | 
   | NEWRING = ATFITCHROM(RING,...,'dp',DP)
   |    Specify off-momentum if radiation is OFF (default 0)
   | 
   | NEWRING = ATFITCHROM(RING,...,'dct',DCT)
   |    Specify path lengthening if radiation is OFF (default 0)
   | 
   |  See also ATFITTUNE

.. py:function:: atrotatelattice

   | ATROTATELATTICE    Circularly shift the lattice elements
   | 
   | NEWRING = ATROTATELATTICE(RING,INDEX)
   |    Return a new lattice such that RING(INDEX) is the first element

.. py:function:: attiltelem

   | ATTILTELEM sets new rotation parameters
   |   NEWELEM = ATTILTELEM(OLDELEM,PSI)
   | 
   |   PSI - rotation angle in radians
   |    POSITIVE ROTS corresponds to a CORKSCREW (right)
   |    rotation of the ELEMENT looking in the direction of the beam.
   |    (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
   |    The rotation matrixes are stored in fields R1 and R2
   |    R1 = [cos(PSI) sin(PSI); -sin(PSI) cos(PSI)]
   |    R2 = R1'
   | 
   |    NEWELEM=ATTILTELEM(OLDELEM,PSI,'RelativeTilt')
   |    the rotation is added to the previous misalignment
   | 
   |   See also ATSETTILT, ATSHIFTELEM, ATMODELEM

.. py:function:: atmakefielderrstruct

   |  MAKERNDFIELDERRS will create a field error data structure
   | 
   |  class is an element class, "Quadrupole", or "Sextupole"

.. py:function:: isatelem

   | ISATELEM tests if an input argument is a valid AT element.
   | 
   |   A valid AT element is a MATLAB structure with required
   |    fields 'Length', 'PassMethod', and a set of data fields,
   |    specific to the PassMethod used.
   | 
   |   [TEST, ERRORSTR] = ISATELEM(ELEM)
   |                    = ISATELEM(ELEM, 'display')
   | 
   |   TEST     - test result,  1 = valid AT element
   |   ERRORSTR - multi-line error message
   | 
   |   See also PASSMETHOD, ATELEM

.. py:function:: atwritepy

   | ATWRITEPY Creates pyAT lattice from a Matlab lattice
   | 
   | PYRING=ATWRITEPY(RING)
   |    return the python lattice object
   | 
   | PYRING=ATWRITEPY(RING,'file',FILENAME)
   |    In addition, store the python lattice object in FILENAME. The ring can
   |    be retrieved in python with the commands:
   | 
   |    >>> with open(<filename>,'rb') as fid:
   |    ...  ring=pickle.load(fid)

.. py:function:: atsplitelem

   | ATSPLITELEM Creates a line by inserting one or more elements into a base element
   | 
   | LINE=ATSPLITELEM(BASEELEM,FRAC1,ELEM1[,FRAC2,ELEM2...])
   |    Each inserted element is associated with a location given by 0<=FRAC<=1
   |    LINE is a cell array containing the sequence of resulting elements
   | 
   |  FRACi may be a scalar or a line vector of locations. ELEMi must be a
   |  single element, or a cell array of elements with the same length as FRACi.
   | 
   |  if FRAC = 0, the element is inserted before BASEELEM (no splitting)
   |  if FRAC = 1, the element is inserted after BASEELEM (no splitting)
   | 
   |  if ELEMi = [], nothing is inserted, only the splitting takes place
   | 
   |  ATSPLITELEM will split BASEELEM in elements with negative lengths if
   |  necessary. Elements with length 0 are not generated. When splitting
   |  dipoles, bending angles are distributed as lengths, and face angles are
   |  set to zero at splitting points.
   | 
   | 
   |  Examples:
   | 
   | >> dr=atdrift('DR',2);       % Insert quadrupoles inside a drift
   | >> qf=atquadrupole('QF',0.1,0.5);
   | >> qd=atquadrupole('QD',0.1,-0.5);
   | >> line=atsplitelem(dr,0.2,qf,0.8,qd);
   | 
   | >> mk=atmarker('MK');
   | >> line=atsplitelem(qf,0,mk);   % Insert a marker before a quadrupole
   | 
   | >> line=atsplitelem(qf,0.5,[]); % Split a quadrupole in two halves
   | 
   |  See also ATINSERTELEMS ATDIVELEM

.. py:function:: splitdrift

   | SPLITDRIFT inserts an element into a drift space
   | 
   |  SPLITDRIFT(DRIFTPOS, SPLIT) inserts a marker (zero-length) element
   |    at distance SPLIT ( 0 < SPLIT < 1) into a drift space
   |    located at DRIFTPOS in THERING
   | 
   |  SPLITDRIFT(DRIFTPOS, SPLIT, ELEMSTRUCCTURE) inserts a marker (zero-length) element
   |    at distance SPLIT ( 0 < SPLIT < 1) into a drift space
   |    located at DRIFTPOS in THERING
   | 
   |  Number of elements in the RING is thus increased by 2
   |  SPLIT (controls the position of the split
   |  L1 = L0*SPLIT
   |  L2 = L0(1-SPLIT)
   |   where L0 is the length of the original DRIFT
   | 
   |  See also: MERGEDRIFT

.. py:function:: atCheckRingProperties

   | ATCHECKRINGPROPERTIES Get the ring properties if existing
   | 
   |  PROPERTIES=ATCHECKRINGPROPERTIES(ring)
   | 	Extract data from the RingParam element of the lattice, if present,
   | 	otherwise return an empty structure
   | 
   |   See also ATGETRINGPROPERTIES, ATSETRINGPROPERTIES

.. py:function:: atshiftelem

   | ATSHIFTELEM set new displacement parameters
   | NEWELEM=ATSHIFTELEM(OLDELEM,DX,DZ)
   | 
   | DX:	horizontal element displacement
   | DY:	Vertical element displacement
   | 
   | 
   | NEWELEM=ATSHIFTELEM(OLDELEM,DX,DY,'RelativeShift')
   | the shift is added to the previous misalignment of the element
   | 
   | See also: atsetshift, attiltelem, atmodelem

.. py:function:: atfastring

   | ATFASTRING Generate simplified AT structures
   | 
   |  The given ring structure is modified so that cavities are moved to the
   |  beginning and the rest of the ring is replaced by a linear 6x6 transfer
   |  matrix followed by a non-linear element providing tune shifts with
   |  amplitude and momentum.
   | 
   |    [FASTRING,FASTRINGRAD]=ATFASTRING(RING)
   | 
   | RING:          original AT structure, with no RF and no radiation.
   | 
   | FASTRING:      Structure containing unchanged cavities moved to the
   |                beginning, a linear 6x6 matrix and a  non-linear element
   |                simulating linear chromaticities and tune shift with
   |                amplitudes
   | 
   | FASTRINGRAD:   Structure containing unchanged cavities moved to the
   |                beginning, a diffusion element, a linear 6x6 transfer
   |                matrix and a non-linear element simulating linear
   |                chromaticities and tune shift with amplitudes
   | 
   |    [FASTRING,FASTRINGRAD]=ATFASTRING(RING,REFPTS)
   | 
   |  The ring is split at the specified locations, and each section is
   |  transformed in the same way as previously described
   | 
   | [FASTRING,FASTRINGRAD]=ATFASTRING(RING,'Plot') plots the tune shifts with amplitude

.. py:function:: entrancefields

   | ENTRANCEFIELDS() Return the list of field names affecting the element entrance
   | 
   |  Optional arguments:
   |  'KeepAxis', if present, rotations translations are excluded from list
   | 
   | see also: exitfields atdivelem

.. py:function:: getcellstruct

   | GETCELLSTRUCT retrieves the field values MATLAB cell array of structures
   | 
   |  VALUES = GETCELLSTRUCT(CELLARRAY,'field',INDEX,M,N)
   | 
   |  VALUES = GETCELLSTRUCT(CELLARRAY,'field',INDEX,M) can be used
   |    for one dimensional vectors
   | 
   |  VALUES = GETCELLSTRUCT(CELLARRAY,'field',INDEX) is the same as
   |    GETCELLSTRUCT(CELLARRAY,'field',index,1,1) if the field data
   |    is a scalar
   | 
   |  VALUES = GETCELLSTRUCT(CELLARRAY,'field',INDEX) is a MATLAB cell array
   |  	 of strings if specified fields contain strings.
   | 
   |  See also ATGETFIELDVALUES SETCELLSTRUCT FINDCELLS

.. py:function:: atloadfielderrs

   |  ATLOADFIELDERRS will load a field error structure into a ring
   |  field error structure has the fields:
   |  elemindex: element indices in ring to impact
   |  Nval: reference mpole #, e.g. 2 for Quad, 3 for Sextupole
   |  nval: multipole index to change
   |  Bval: Value of relative normal coefficient
   |  Aval: Value of relative skew coefficient
   |  radius: reference radius used to compute error

