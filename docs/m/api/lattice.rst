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
     - CONVERTERS
   * - :ref:`paramgroup_module`
     - PARAMGROUP
   * - :ref:`element_creation_module`
     - ELEMENT_CREATION
   * - :ref:`survey_module`
     - SURVEY

.. rubric:: Classes


.. list-table::

   * - :class:`atparticle`
     - 

.. rubric:: Functions


.. list-table::

   * - :func:`at2py`
     - ELSTR=(ELEM) convert AT element tp pyat
   * - :func:`at2str`
     - Makes the string representation of an AT element
   * - :func:`atCheckRingProperties`
     - Get the ring properties if existing
   * - :func:`atGetRingProperties`
     - Get the ring properties
   * - :func:`atSetRingProperties`
     - atSetRingProperties	Add or modify properties of the lattice
   * - :func:`ataddmpolecomppoly`
     - ataddmpolecomppoly adds a multipole component to an existing polynomial,
   * - :func:`ataddmpoleerrors`
     - ataddrandmpole adds a random multipole component to all elements of type
   * - :func:`atdivelem`
     - Divide an element into pieces
   * - :func:`atelem`
     - makes a new AT element structure from another element,
   * - :func:`atfastring`
     - Generate simplified AT structures
   * - :func:`atfitchrom`
     - Fit chromaticites by scaling 2 sextupole families
   * - :func:`atfittune`
     - Fit linear tunes by scaling 2 quadrupole families
   * - :func:`atgetcells`
     - performs a search on MATLAB cell arrays of structures
   * - :func:`atgetfieldvalues`
     - retrieves the field values AT cell array of elements
   * - :func:`atguessclass`
     - Tries to determine the class of an element
   * - :func:`atindex`
     - extracts the information about element families and
   * - :func:`atinsertelems`
     - Insert elements at given locations in a line
   * - :func:`atloadfielderrs`
     - will load a field error structure into a ring
   * - :func:`atloadlattice`
     - load a lattice from a file
   * - :func:`atmaincavities`
     - Get the fundamental mode cavities
   * - :func:`atmakefielderrstruct`
     - MAKERNDFIELDERRS will create a field error data structure
   * - :func:`atreduce`
     - Remove useless elements from an AT structure
   * - :func:`atrotatelattice`
     - Circularly shift the lattice elements
   * - :func:`atsbreak`
     - Insert markers at given s positions in a lattice
   * - :func:`atsetRFCavity`
     - Set the RF Cavity with the passmethod RFCavityPass.
   * - :func:`atsetfieldvalues`
     - sets the field values of MATLAB cell array of structures
   * - :func:`atsetshift`
     - sets the misalignment vectors
   * - :func:`atsettilt`
     - sets the entrance and exit rotation matrices
   * - :func:`atshiftelem`
     - set new displacement parameters
   * - :func:`atsimplering`
     - 
   * - :func:`atsplitelem`
     - Creates a line by inserting one or more elements into a base element
   * - :func:`attiltelem`
     - sets new rotation parameters
   * - :func:`atwritejson`
     - Create a JSON file to store an AT lattice
   * - :func:`atwritem`
     - Creates a .m file to store an AT structure
   * - :func:`atwritepy`
     - Creates pyAT lattice from a Matlab lattice
   * - :func:`buildlat`
     - places elements from FAMLIST into cell array THERING
   * - :func:`combinebypassmethod`
     - combines adjacent elements that have the same specified pass method
   * - :func:`combinelinear45`
     - combines adjacent  elements that use 4-by-5 PassMethods
   * - :func:`entrancefields`
     - () Return the list of field names affecting the element entrance
   * - :func:`exitfields`
     - () Return the list of field names affecting the element exit
   * - :func:`findcells`
     - performs a search on MATLAB cell arrays of structures
   * - :func:`findtags`
     - looks for string matches in 'Tag' field of AT lattice elements
   * - :func:`getcellstruct`
     - retrieves the field values MATLAB cell array of structures
   * - :func:`insertelem0`
     - - quick and dirty:
   * - :func:`insertindrift`
     - inserts one or more elements into a drift element
   * - :func:`isatelem`
     - tests if an input argument is a valid AT element.
   * - :func:`isatlattice`
     - tests if an input argument is a valid AT lattice.
   * - :func:`mergedrift`
     - removes a lattice element and merges the two adjacent drift spaces
   * - :func:`mvelem`
     - Move an element
   * - :func:`mvfield`
     - Move fields from one structure to another
   * - :func:`rmelem0`
     - removes elements of length 0 from the accelerator lattice
   * - :func:`setcellstruct`
     - sets the field values of MATLAB cell array of structures
   * - :func:`setshift`
     - sets the misalignment vectors T1, T2 for elements
   * - :func:`settags`
     - sets the 'Tag' field in AT lattice elements
   * - :func:`settilt`
     - sets the entrance and exit misalignment matrixes
   * - :func:`splitdrift`
     - inserts an element into a drift space
   * - :func:`symplectify`
     - symplectify makes a matrix more symplectic

.. py:function:: atparticle


.. py:function:: at2py

   | ELSTR=(ELEM) convert AT element tp pyat
   
   | **at2py** Creates a pyat element

.. py:function:: at2str(elem))

   | Makes the string representation of an AT element
   
   | **at2str** Creates a string such that EVAL(**at2str(elem))** recreates an
   | identical element.
   
   |   INPUTS
   |   1. elem - Elem to write
   
   |   OUTPUTS
   |   1. elsstr - String given the AT constructor of the element
   
   | See also :func:`atwritem`, :func:`atwritepy`

.. py:function:: atCheckRingProperties(ring)

   | Get the ring properties if existing
   
   |  **properties=atCheckRingProperties(ring)**
   | 	Extract data from the RingParam element of the lattice, if present,
   | 	otherwise return an empty structure
   
   | See also :func:`atgetringproperties`, :func:`atsetringproperties`

.. py:function:: atGetRingProperties(ring,'param1','param2',...)

   | Get the ring properties
   
   |  **[v1,v2,...]=atGetRingProperties(ring,'param1','param2',...)**
   |    Extract lattice properties. Extract from the RingParam element of the
   |    lattice if present, or from the lattice elements. The number of outputs
   |    corresponds to the number of property names in input.
   
   |  RING:             Ring structure
   
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
   
   |  Properties for the full ring ('periods' x cells):
   |    'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
   |    'harmonic_number'          "       "
   |    'Circumference'         Ring circumference [m] (cell_length * Periodicity)
   |    'rf_voltage'            RF voltage [V] (cell_rf_voltage * Periodicity)
   |    'revolution_frequency'  Revolution frequency [Hz] (cell_revolution_frequency / Periodicity)
   
   |  Properties for one cell:
   |    'cell_harmnumber'       Harmonic number (for 1 cell)
   |    'cell_length'           Cell length [m]
   |    'cell_rf_voltage'       RF voltage per cell [V] (main cavities)
   |    'cell_revolution_frequency' Revolution frequency [Hz] (for 1 cell)
   
   |  Custom properties may be added with atSetRingProperties. Custom property
   |  names are case dependent.
   
   |  Example:
   |  >> **[energy, beta] = atGetRingProperties(ring,'energy','beta')**;
   
   |  **properties=atGetRingProperties(ring)**
   
   |  RING:             Ring structure
   |  PROPERTIES:       Structure with fields:
   |    'FamName'               Name of the lattice
   |    'Energy'                Ring energy [eV]
   |    'Periodicity'           Number of periods to build a full ring
   |    'Particle'              particle (Particle object)
   |    'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
   |    'cell_harmnumber'       Harmonic number (for 1 cell)
   |    'cavpts'                Location of the main cavities
   
   |  **properties=atGetRingProperties(ring,'all')**
   |    The PROPERTIES structure contains all the properties described above,
   |    both standard and custom.
   
   |  For fast access, some ring properties are stored in a RingParam element
   |  ideally located in the 1st position of the lattice. Without such element,
   |  the properties are deduced from the lattice contents. This is much slower
   |  and **atGetRingProperties** displays a warning indicating how to add the
   |  RingParam element:
   
   | >>ring=atSetRingProperties(ring)
   
   | See also :func:`atsetringproperties`

.. py:function:: atSetRingProperties(ring [,key,value]...)

   | atSetRingProperties	Add or modify properties of the lattice
   
   |  **newring=atSetRingProperties(ring [,key,value]...)**
   |    Add or modify the attributes of the RingParam element of the lattice,
   |    Insert a new RingParam element if necessary
   
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
   
   |  Properties for the full ring ('periods' x cells):
   |    'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
   |    'harmonic_number'          "       "
   |    'rf_voltage'            RF voltage [V] (cell_rf_voltage * Periodicity)
   
   |  Properties for one cell:
   |    'cell_harmnumber'       Harmonic number (for 1 cell)
   |    'cell_rf_voltage'       RF voltage per cell [V] (main cavities)
   
   |  Additional custom fields may be added. They can be retrieved by
   |  atGetRingProperties and are saved in files. Custom field names are
   |  case-dependent.
   
   |  For fast access, some ring properties are stored in a RingParam element
   |  ideally located in the 1st position of the lattice. If there is no such
   |  element, **atSetRingProperties** will add it.
   
   | See also :func:`atgetringproperties`

.. py:function:: ataddmpolecomppoly(polynom,index,strength,radius)

   | ataddmpolecomppoly adds a multipole component to an existing polynomial,
   | scaling it by the reference multipole value and a radius to set the length
   | scale
   | **[polynomout] = ataddmpolecomppoly(polynom,index,strength,radius)**
   
   | Polynom = given field polynomial that we want to add a component to
   | refindex = index of magnet: 1 for dipole, 2 for quadrupole
   | newindex: index of Multipole to add
   | strength: strength of the multipole component at the given radius
   | radius: reference radius for defining the absolute strength
   |  B^(N)_(n) = radius^(n-N)*b_n/b_N
   
   | See also :func:`attiltelem`, :func:`atshiftelem`, :func:`atmodelem`

.. py:function:: ataddmpoleerrors

   | ataddrandmpole adds a random multipole component to all elements of type
   | 'type' where type can be 'dipole', 'quadrupole', or 'sextupole'
   
   | [newring] = ATRANDMPOLE(ring,type,newindex,strength,radius)
   
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
   
   | See also :func:`ataddmpolecomppoly`, :func:`attiltelem`, :func:`atshiftelem`

.. py:function:: atdivelem(elem,frac)

   | Divide an element into pieces
   
   | **line=atdivelem(elem,frac)**
   
   | ELEM:  Element to be divided
   | FRAC:  Length of each piece, as a fraction of the total length
   
   | LINE:  Cell array
   
   |  The sum of frac elements does not need to be 1, however for bending
   |  magnets, the length will be divided according to FRAC, but the total
   |  bending angle will be divided according to FRAC/SUM(FRAC) so that the
   |  total bending angle is kept.
   
   |  Example:
   
   | >> qf=atquadrupole('QF',0.1,0.5);
   | >> **line=atdivelem(qf,[0.5;0.5])**; % Split a quadrupole in two halves
   
   |  Optional arguments:
   |  'KeepAxis', if present, rotations translations are kept at all slices
   
   | See also :func:`atinsertelems`, :func:`atslice`, :func:`atsplitelem`, :func:`entrancefields`, :func:`exitfields`

.. py:function:: atelem(elem,'field1','value1','field2', 'value2', ...)

   | makes a new AT element structure from another element,
   |  standard AT type, or a template on file. Existing fields are overridden
   |  with new values. New fields are added and initialized.
   
   |  **newelem = atelem(elem,'field1','value1','field2', 'value2', ...)**
   |   ELEM can be 1) another element structure
   |               2) One of the standard AT types
   |                  'drift'
   |                  'quadrupole'
   |                  'sextupole'
   |                  'bend','rbend','sbend'
   |                  'marker'
   |                  'corrector'
   |                   ...

.. py:function:: atfastring(ring)

   | Generate simplified AT structures
   
   |  The given ring structure is modified so that cavities are moved to the
   |  beginning and the rest of the ring is replaced by a linear 6x6 transfer
   |  matrix followed by a non-linear element providing tune shifts with
   |  amplitude and momentum.
   
   |    **[fastring,fastringrad]=atfastring(ring)**
   
   | RING:          original AT structure, with no RF and no radiation.
   
   | FASTRING:      Structure containing unchanged cavities moved to the
   |                beginning, a linear 6x6 matrix and a  non-linear element
   |                simulating linear chromaticities and tune shift with
   |                amplitudes
   
   | FASTRINGRAD:   Structure containing unchanged cavities moved to the
   |                beginning, a diffusion element, a linear 6x6 transfer
   |                matrix and a non-linear element simulating linear
   |                chromaticities and tune shift with amplitudes
   
   |    **[fastring,fastringrad]=atfastring(ring,refpts)**
   
   |  The ring is split at the specified locations, and each section is
   |  transformed in the same way as previously described
   
   | **[fastring,fastringrad]=atfastring(ring,'plot')** plots the tune shifts with amplitude

.. py:function:: atfitchrom(ring,newchrom,sextfamily1,sextfamily2)

   | Fit chromaticites by scaling 2 sextupole families
   
   |  **newring = atfitchrom(ring,newchrom,sextfamily1,sextfamily2)**
   |  **newring = atfitchrom(ring,dpp,newchrom,sextfamily1,sextfamily2)**
   
   | RING:          Cell array
   | DPP:           Optional momentum deviation (default 0)
   | NEWCHROM:      Desired non-normalized chromaticities
   | SEXTFAMILY1:   1st sextupole family
   | SEXTFAMILY2:   2nd sextupole family
   
   | SEXTFAMILY may be:
   |    string: Family name
   |    logical array: mask of selected elements in RING
   |    Numeric array: list of selected elements in RING
   |    Cell array: All elements selected by each cell
   
   | **newring = atfitchrom(ring,...,'dpstep',dpstep)**
   |    dpstep is the momentum step applied to compute the chromaticity.
   |    The default is the DPStep global option, which defaults to 3.0e-6
   
   | **newring = atfitchrom(ring,...,'hstep',hstep)**
   |    hstep is the sextupole strength applied to build the jacobian [m^-3].
   |    Default: 1.0e-5
   
   | **newring = atfitchrom(ring,...,'dp',dp)**
   |    Specify off-momentum if radiation is OFF (default 0)
   
   | **newring = atfitchrom(ring,...,'dct',dct)**
   |    Specify path lengthening if radiation is OFF (default 0)
   
   | See also :func:`atfittune`

.. py:function:: atfittune(ring,newtunes,quadfamily1,quadfamily2)

   | Fit linear tunes by scaling 2 quadrupole families
   
   |  **newring = atfittune(ring,newtunes,quadfamily1,quadfamily2)**
   |  **newring = atfittune(ring,dpp,newtunes,quadfamily1,quadfamily2)**
   
   | RING:          Cell array
   | DPP:           Optional momentum deviation (default 0)
   | NEWTUNES:      Desired tune values
   | QUADFAMILY1:   1st quadrupole family
   | QUADFAMILY2:   2nd quadrupole family
   
   | QUADFAMILY may be:
   |    string: Family name
   |    logical array: mask of selected elements in RING
   |    Numeric array: list of selected elements in RING
   |    Cell array: All elements selected by each cell
   
   | **newring = atfittune(ring,...,'useintegerpart')** With this flag, the
   |    function fits the tunes to the total values of NEWTUNES, including
   |    the integer part.
   |    With this option the function is substantially slower!
   
   | **newring = atfittune(ring,...,'kstep',kstep)**
   |    kstep is the quadrupole strength applied to build the jacobian [m^-2].
   |    Default: 1.0e-6
   
   | **newring = atfittune(ring,...,'dp',dp)**
   |    Specify off-momentum if radiation is OFF (default 0)
   
   | **newring = atfittune(ring,...,'dct',dct)**
   |    Specify path lengthening if radiation is OFF (default 0)
   
   | See also :func:`atfitchrom`

.. py:function:: atgetcells(ring, 'field')

   | performs a search on MATLAB cell arrays of structures
   
   |  **ok = atgetcells(ring, 'field')**
   |    returns indexes of elements that have a field named 'field'
   
   |  **ok = atgetcells(ring, 'field', value1...)**
   |    returns indexes of elements whose field 'field'
   |    is equal to VALUE1, VALUE2, ... or VALUEN. Where VALUE can either be
   |    character strings or a number. If its a character string REGULAR
   |    expressions can be used.
   
   |  **ok = atgetcells(ring, 'field', @testfunction, args...)**
   |    Uses the user-defined TESTFUNCTION to select array elements
   |    TESTFUNCTION must be of the form:
   |        OK=TESTFUNTION(ATELEM,FIELDVALUE,ARGS...)
   
   |  OK is a logical array with the same size as RING, refering to matching
   |  elements in RING
   
   | See also :func:`atgetfieldvalues`, :func:`atsetfieldvalues`, :func:`findcells`, :func:`regexpi`

.. py:function:: atgetfieldvalues(ring,'field')

   | retrieves the field values AT cell array of elements
   
   |  **values = atgetfieldvalues(ring,'field')** extracts the values of
   |  the field 'field' in all the elements of RING
   
   |  **values = atgetfieldvalues(ring,index,'field')** extracts the values of
   |  the field 'field' in the elements of RING selected by INDEX
   
   |  if RING{I}.FIELD is a numeric scalar
   |     VALUES is a length(INDEX) x 1 array
   |  otherwise
   |     VALUES is a length(INDEX) x 1 cell array
   
   | **values = atgetfieldvalues(...,'default',default_value)** Uses default_values
   |    if the required field is not existing
   
   |  More generally **atgetfieldvalues(ring,index,subs1,subs2,...)** will call
   |   GETFIELD(RING{I},subs1,subs2,...) for I in INDEX
   
   |  Examples:
   
   |  **v=atgetfieldvalues(ring,1:10,'polynomb')** is a 10x1  cell array
   |  such that V{I}=RING{I}.PolynomB for I=1:10
   
   |  **v=atgetfieldvalues(ring(1:10),'polynomb',{1,2})** is a 10x1 array
   |  such that V(I)=RING{I},PolynomB(1,2)
   
   
   | See also :func:`atsetfieldvalues`, :func:`atgetcells`, :func:`getcellstruct`, :func:`findcells`

.. py:function:: atguessclass(atelem)

   | Tries to determine the class of an element
   | **atclass=atguessclass(atelem)** Tries to determine the class of an element
   
   |   INPUTS
   |   1. elem - AT element
   
   
   |   NOTES
   |   1. **atclass=atguessclass(atelem,'useclass')**
   |   By default, **atguessclass** will default "useless" elements (PolynopmB==0)
   |   to 'Drift' or 'Marker', depending on 'Length'. When specifying
   |   'UseClass', it it will preserve the 'Class' field for those elements.
   
   | See also :func:`atwritem`

.. py:function:: atindex(lattice)

   | extracts the information about element families and
   |  indexes from AT lattice
   
   | **ati=atindex(lattice)**
   |   returns a srtucture with fields named after element family
   |   containing an array of their AT indexes;
   
   |    ATI.QF = [...]
   |    ATI.QD = [...];
   |    ...
   
   | **ati=atindex()** Uses the global variable THERING
   
   | See also :func:`atgetcells`

.. py:function:: atinsertelems(line,refpts,frac1,elem1[,frac2,elem2...])

   | Insert elements at given locations in a line
   
   | **newline=atinsertelems(line,refpts,frac1,elem1[,frac2,elem2...])**
   |    a new line is created by inserting elements at each location specified
   |    by REFPTS.
   
   | LINE:      Cell array of structures
   | REFPTS:    Insertion points (index list or logical mask)
   | FRAC:      Location of the inserted element ELEM within LINE{REFPTS(i)}
   |            0<=FRAC<=1
   |  if FRAC = 0, ELEM is inserted before LINE{REFPTS(i)} (no splitting)
   |  if FRAC = 1, ELEM is inserted after LINE{REFPTS(i)} (no splitting)
   |  if FRAC = NaN, LINE{REFPTS(i)} is replaced by ELEM (no check for identical length)
   |  if ELEM = [], nothing is inserted, only the splitting takes place
   
   |  FRAC and ELEM must be scalars or array of the same size as REFPTS
   
   | See also :func:`atsplitelem`, :func:`atdivelem`

.. py:function:: atloadfielderrs

   | will load a field error structure into a ring
   |  field error structure has the fields:
   |  elemindex: element indices in ring to impact
   |  Nval: reference mpole #, e.g. 2 for Quad, 3 for Sextupole
   |  nval: multipole index to change
   |  Bval: Value of relative normal coefficient
   |  Aval: Value of relative skew coefficient
   |  radius: reference radius used to compute error

.. py:function:: atloadlattice(filepath)

   | load a lattice from a file
   
   | **lattice=atloadlattice(filepath)**
   |    Create an AT lattice (cell array of AT elements) from FILEPATH
   
   | **lattice=atloadlattice(filepath,'matkey',variable_name,...)**
   |    Give the name of a variable containing the lattice in .mat files
   |    containing several variables. By default, if the .mat file contains a
   |    single variable, it will be loaded. Otherwise **atloadlattice** will look
   |    in order for 'ring', 'lattice' and 'RING'.
   
   | **lattice=atloadlattice(filepath,key,value...)**
   |    (Key,value) pairs are added to the RingProperties of the lattice or
   |    overload existing ones.
   |    Standard keys include: FamName, Energy, Periodicity, HarmNumber, Particle
   |    Custom properties are allowed
   
   | Allowed file types:
   |    .mat    Binary Matlab file. If the file contains several variables,
   |            the variable name must be specified using the 'matkey' keyword.
   
   |    .m      Matlab function. The function must output a valid AT structure.
   |    .json   JSON file
   
   | see also atwritem, atwritejson

.. py:function:: atmaincavities(ring)

   | Get the fundamental mode cavities
   
   | **idcav=atmaincavities(ring)**
   |    Return a refpts pointing to fundamental mode cavities
   |    (cavities with the lowest frequency)

.. py:function:: atmakefielderrstruct

   | MAKERNDFIELDERRS will create a field error data structure
   
   |  class is an element class, "Quadrupole", or "Sextupole"

.. py:function:: atreduce(ring)

   | Remove useless elements from an AT structure
   | **newring=atreduce(ring)**
   
   |  Remove elements with PassMethod='IdentityPass' and merges adjacent
   |  similar elements
   
   | **newring=atreduce(ring,refpts)**
   | 	When merging similar elements, keep REFPTS intact.
   
   | **[newring,kept]=atreduce(...)**
   | 	Returns the index of kept elements so that NEWRING=OLDRING(KEPT)
   
   | **[newring,kept,newrefpts]=atreduce(ring,refpts)**
   | 	Returns in addition the updated list of reference points.
   

.. py:function:: atrotatelattice(ring,index)

   | Circularly shift the lattice elements
   
   | **newring = atrotatelattice(ring,index)**
   |    Return a new lattice such that RING(INDEX) is the first element

.. py:function:: atsbreak(ring,spos)

   | Insert markers at given s positions in a lattice
   
   | **[newring,refpts]=atsbreak(ring,spos)** build a new lattice with breakpoints
   
   | RING:  AT structure
   | SPOS:  Poition of the desired breakpoints
   
   | NEWRING:   Modified AT structure
   | REFPTS:    Index of breakpoints
   
   | See also :func:`atinsertelems`, :func:`atsplitelem`, :func:`atdivelem`

.. py:function:: atsetRFCavity(ring, rfv, radflag, harmnumber, deltafreq)

   | Set the RF Cavity with the passmethod RFCavityPass.
   |   RFCavityPass allows to change the energy of the beam changing the
   |   frequency of the cavity.
   
   |    **newring = atsetRFCavity(ring, rfv, radflag, harmnumber, deltafreq)**
   |    sets the synchronous phase of the cavity, the voltage, the harmonic
   |    number and the PassMethod.
   
   |   INPUTS
   |   1. ring       - Ring structure
   |   2. rfv        - RF-voltage in volts
   |   3. radflag    - 0/1: activat/desactivate radiation (atradon/atradoff)
   |   4. HarmNumber - Harmonic number
   |   5. DeltaFreq  - Frequency shift in Hz
   
   |   OUTPUTS
   |   1. newring - update ring structure with nw RF parameters
   
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
   
   |   EXAMPLES
   
   |    1. normal use:
   |    **newring = atsetRFCavity( ring, 6e6, 1, 992, 0 )**
   
   |    2. for off-energy simulations:
   |    **newring = atsetRFCavity( ring, 6e6, 1, 992, 100 )**
   
   | See also :func:`atsetcavity`, :func:`atradon`, :func:`atradoff`, :func:`atgetu0`

.. py:function:: atsetfieldvalues(ring,...)

   | sets the field values of MATLAB cell array of structures
   
   |    Note that the calling syntax must be in the form of assignment:
   |    **ring = atsetfieldvalues(ring,...)**
   |    MATLAB does not modify variables that only appear on the right
   |    hand side as arguments.
   
   | **newring=atsetfieldvalues(ring,'field',values)**
   |    In this mode, the function will set values on all the elements of RING
   
   | **newring=atsetfieldvalues(ring,index,'field',values)**
   |    In this mode, the function will set values on the elements of RING
   |    specified by INDEX, given as a list of indices or as a logical mask
   
   | **newring=atsetfieldvalues(ring,'field',valuestruct)**
   |    In this mode, the function will set values on the elements of RING
   |    whose family names are given by the field names of VALUESTRUCT
   
   | **newring=atsetfieldvalues(ring,ringindex,....,valuestruct)**
   |    As in the previous mode, the function will set values on the elements
   |    of RING whose family names are given by the field names of VALUESTRUCT.
   |    But RINGINDEX=atindex(RING) is provided to avoid multiple computations.
   
   |  Field selection
   |  ---------------------------------------------------------
   |  NEWRING = ATSETFIELD(RING,'field',VALUES)
   |    For each I=1:length(RING), set RING{I}.FIELD=value
   
   |  NEWRING = ATSETFIELD(RING,'field',{M,N},VALUES)
   |    For each I=1:length(RING), set RING{I}.FIELD(M,N)=value
   
   |  More generally,
   |  NEWRING = ATSETFIELD(RING,subs1,subs2,...,VALUES)
   |    For each I=1:length(RING), SETFIELD(RING{I},subs1,subs2,...,value)
   
   |  The first dimension of VALUES must be either length(INDEX) or 1 (the value
   |  will be repeated for each element). For a vector to be repeated, enclose
   |  it in a cell array.
   
   |  Value format
   |  ---------------------------------------------------------
   |  Cell array VALUES
   |  -----------------
   |  Mx1 or 1xM cell array : one cell per element
   |  1x1 cell array : cell 1 is affected to all selected elements
   
   |  Character array VALUES
   |  ---------------------
   |  1xN char array (string) : the string as affected to all selected elements
   |  MxN char array : one row per element
   
   |  Numeric array VALUES
   |  --------------------
   |  1x1 (scalar) : the value is affected to all selected elements
   |  Mx1 (column) : one value per element
   |  MxN (matrix) : one row affected per element. If M==1, the single row
   |                 is affected to all selected elements
   |  To assign column vectors to parameters, use a cell array VALUES, with one
   |  column per cell
   
   | See also :func:`atgetfieldvalues`, :func:`atgetcells`, :func:`setcellstruct`, :func:`findcells`

.. py:function:: atsetshift(ring, elemindex, dx, dy)

   | sets the misalignment vectors
   
   |  **ring=atsetshift(ring, elemindex, dx, dy)** sets the entrance and exit misalignment vectors
   |   of one element or a group of elements in the globally defined lattice THERING.
   
   |  ELEMINDEX contains indexes of elements to be misaligned
   |  DX, DY are displacements of the ELEMENT
   |   so the coordinate transformation on the particle at entrance is
   | 	X  ->  X-DX
   |    Y  ->  Y-DY
   
   |  **ring=atsetshift(ring, elemindex, dx, dy, 'relativeshift')**
   |  adds the shifts DX and DY to the the previous misalignments
   |  of the elements
   
   |  **atsetshift(elemindex, dx, dy)** Uses the global variable THERING
   
   | See also :func:`atsettilt`

.. py:function:: atsettilt(ring,elemindex, psi)

   | sets the entrance and exit rotation matrices
   |  of an element or a group of elements in THERING
   
   |  **ring=atsettilt(ring,elemindex, psi)**
   |  ELEMINDEX contains indexes of elements to be rotated
   |  PSI - angle(s) of rotation in RADIANS
   |    POSITIVE PSI corresponds to a CORKSCREW (right)
   |    rotation of the ELEMENT looking in the direction of the beam.
   |    (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
   |    The misalgnment matrixes are stored in fields R1 and R2
   |    R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
   |    R2 = R1'
   
   |  **ring=atsettilt(ring,elemindex,psi,'relativetilt')**
   |  the rotation is added to the previous misalignment
   
   |  **atsettilt(elemindex, psi)** Uses the global variable THERING
   
   | See also :func:`atsetshift`

.. py:function:: atshiftelem(oldelem,dx,dz)

   | set new displacement parameters
   | **newelem=atshiftelem(oldelem,dx,dz)**
   
   | DX:	horizontal element displacement
   | DY:	Vertical element displacement
   
   
   | **newelem=atshiftelem(oldelem,dx,dy,'relativeshift')**
   | the shift is added to the previous misalignment of the element
   
   | See also :func:`atsetshift`, :func:`attiltelem`, :func:`atmodelem`

.. py:function:: atsimplering


.. py:function:: atsplitelem(baseelem,frac1,elem1[,frac2,elem2...])

   | Creates a line by inserting one or more elements into a base element
   
   | **line=atsplitelem(baseelem,frac1,elem1[,frac2,elem2...])**
   |    Each inserted element is associated with a location given by 0<=FRAC<=1
   |    LINE is a cell array containing the sequence of resulting elements
   
   |  FRACi may be a scalar or a line vector of locations. ELEMi must be a
   |  single element, or a cell array of elements with the same length as FRACi.
   
   |  if FRAC = 0, the element is inserted before BASEELEM (no splitting)
   |  if FRAC = 1, the element is inserted after BASEELEM (no splitting)
   
   |  if ELEMi = [], nothing is inserted, only the splitting takes place
   
   |  **atsplitelem** will split BASEELEM in elements with negative lengths if
   |  necessary. Elements with length 0 are not generated. When splitting
   |  dipoles, bending angles are distributed as lengths, and face angles are
   |  set to zero at splitting points.
   
   
   |  Examples:
   
   | >> dr=atdrift('DR',2);       % Insert quadrupoles inside a drift
   | >> qf=atquadrupole('QF',0.1,0.5);
   | >> qd=atquadrupole('QD',0.1,-0.5);
   | >> **line=atsplitelem(dr,0.2,qf,0.8,qd)**;
   
   | >> mk=atmarker('MK');
   | >> **line=atsplitelem(qf,0,mk)**;   % Insert a marker before a quadrupole
   
   | >> **line=atsplitelem(qf,0.5,[])**; % Split a quadrupole in two halves
   
   | See also :func:`atinsertelems`, :func:`atdivelem`

.. py:function:: attiltelem(oldelem,psi)

   | sets new rotation parameters
   |   **newelem = attiltelem(oldelem,psi)**
   
   |   PSI - rotation angle in radians
   |    POSITIVE ROTS corresponds to a CORKSCREW (right)
   |    rotation of the ELEMENT looking in the direction of the beam.
   |    (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
   |    The rotation matrixes are stored in fields R1 and R2
   |    R1 = [cos(PSI) sin(PSI); -sin(PSI) cos(PSI)]
   |    R2 = R1'
   
   |    **newelem=attiltelem(oldelem,psi,'relativetilt')**
   |    the rotation is added to the previous misalignment
   
   | See also :func:`atsettilt`, :func:`atshiftelem`, :func:`atmodelem`

.. py:function:: atwritejson(ring)

   | Create a JSON file to store an AT lattice
   
   | **js=atwritejson(ring)**
   |    Return the JSON representation of RING as a character array
   
   | **atwritejson(ring, filename)**
   |    Write the JSON representation of RING to the file FILENAME
   
   | **atwritejson(ring, ..., 'compact', true)**
   |    If compact is true, write a compact JSON file (no linefeeds)
   
   | see also atloadlattice

.. py:function:: atwritem(ring)

   | Creates a .m file to store an AT structure
   
   |   INPUTS
   |   1. ring    - Lattice structure (.mat file)
   |   2. filname - output file where to write the lattice as a line by line file using
   |   element constructors
   
   |   EXAMPLES
   |   1. **atwritem(ring)** % Prints the result in the command window
   |   2. **atwritem(ring,filename)** % Prints the result in a file
   
   | See also :func:`at2str`

.. py:function:: atwritepy(ring)

   | Creates pyAT lattice from a Matlab lattice
   
   | **pyring=atwritepy(ring)**
   |    return the python lattice object
   
   | **pyring=atwritepy(ring,'file',filename)**
   |    In addition, store the python lattice object in FILENAME. The ring can
   |    be retrieved in python with the commands:
   
   |    >>> with open(<filename>,'rb') as fid:
   |    ...  ring=pickle.load(fid)

.. py:function:: buildlat

   | places elements from FAMLIST into cell array THERING
   |  in the order given by ineger arry ELIST
   |  to be use in Accelerator Toolbox lattice definition files

.. py:function:: combinebypassmethod(lattice,method,keepindex,refindex)

   | combines adjacent elements that have the same specified pass method
   |  **[newlattice, shiftedkeepindex, shiftedref] = combinebypassmethod(lattice,method,keepindex,refindex)**

.. py:function:: combinelinear45(lattice,keepindex,refindex)

   | combines adjacent  elements that use 4-by-5 PassMethods
   |  **[newlattice, shiftedkeepindex, shiftedref] = combinelinear45(lattice,keepindex,refindex)**

.. py:function:: entrancefields

   | () Return the list of field names affecting the element entrance
   
   |  Optional arguments:
   |  'KeepAxis', if present, rotations translations are excluded from list
   
   | see also: exitfields atdivelem

.. py:function:: exitfields

   | () Return the list of field names affecting the element exit
   
   |  Optional arguments:
   |  'KeepAxis', if present, rotations translations are excluded from list
   
   | see also: entrancefields atdivelem

.. py:function:: findcells(cellarray, 'field')

   | performs a search on MATLAB cell arrays of structures
   
   |  **index = findcells(cellarray, 'field')**
   |    returns indexes of elements that have a field named 'field'
   
   |  **index = findcells(cellarray, 'field', value)**
   |    returns indexes of elements whose field 'field'
   |    is equal to VALUE1, VALUE2, ... or VALUEN. Where VALUE can either be
   |    character strings or a number. If its a character string REGULAR
   |    expressions can be used.
   
   |  Example:
   |    **findcells(thering,'length',0, 0.2)**;  % will match elements of
   |                                           lengths 0 and 0.2
   |    **findcells(thering,'famname','sfa','sda')**;
   
   | See also :func:`atgetcells`, :func:`getcellstruct`, :func:`setcellstruct`, :func:`regexpi`

.. py:function:: findtags(cellarray, matchstr)

   | looks for string matches in 'Tag' field of AT lattice elements
   
   |  **index = findtags(cellarray, matchstr)**
   |    returns indexes of elements that have a field 'Tag'
   |    whose value is a string exactly matching MATCHSTR
   |    or a cell array of strings with one element matching MATCHSTR
   
   | See also :func:`findcells`, :func:`settags`

.. py:function:: getcellstruct(cellarray,'field',index,m,n)

   | retrieves the field values MATLAB cell array of structures
   
   |  **values = getcellstruct(cellarray,'field',index,m,n)**
   
   |  **values = getcellstruct(cellarray,'field',index,m)** can be used
   |    for one dimensional vectors
   
   |  **values = getcellstruct(cellarray,'field',index)** is the same as
   |    **getcellstruct(cellarray,'field',index,1,1)** if the field data
   |    is a scalar
   
   |  **values = getcellstruct(cellarray,'field',index)** is a MATLAB cell array
   |  	 of strings if specified fields contain strings.
   
   | See also :func:`atgetfieldvalues`, :func:`setcellstruct`, :func:`findcells`

.. py:function:: insertelem0(lattice, driftindex, splitlength, elemdata)

   | - quick and dirty:
   |   inserts element(s) of zero length into a drift space
   
   |  **newlattice = insertelem0(lattice, driftindex, splitlength, elemdata)**

.. py:function:: insertindrift(d, qf, 0.5, qd, 2, qf, 3.5)

   | inserts one or more elements into a drift element
   |   and returns a sequence (cell array) of elements  ready to to be used
   |   in AT lattice
   
   |  ELEMSEQ = INSERTELEM(DRIFT0, ELEM1, POS1, ... ELEMN, POSN)
   
   |  EXAMPLE: FODO cell
   
   |  --- 1. Declare elements
   
   |  D  = atelem('drift','Length',4.5);
   |  QF = atelem('quad','Length', 1, 'K',  1.234);
   |  QD = atelem('quad','Length', 1, 'K', -2.345);
   
   |  --- 2. Insert quadrupoles in the drift;
   
   |  **fodocell = insertindrift(d, qf, 0.5, qd, 2, qf, 3.5)**;
   
   | See also :func:`splitelem`

.. py:function:: isatelem(elem)

   | tests if an input argument is a valid AT element.
   
   |   A valid AT element is a MATLAB structure with required
   |    fields 'Length', 'PassMethod', and a set of data fields,
   |    specific to the PassMethod used.
   
   |   **[test, errorstr] = isatelem(elem)**
   |                    = **isatelem(elem, 'display')**
   
   |   TEST     - test result,  1 = valid AT element
   |   ERRORSTR - multi-line error message
   
   | See also :func:`passmethod`, :func:`atelem`

.. py:function:: isatlattice(elem, ['option1',...])

   | tests if an input argument is a valid AT lattice.
   
   |   A valid AT lattice is a MATLAB cell array of valid AT elements
   
   |   **[test, badlist, errorstr] = isatlattice(elem, ['option1',...])**
   
   |   Allowed otions:
   |    'display' - print problems to command window;
   |    'stop'    - return after the first problem is found
   
   |   TEST     - test result,  1 = valid AT element
   |   ERRORSTR - multi-line error message
   
   | See also :func:`isatelem`, :func:`atelem`, :func:`atlattice`

.. py:function:: mergedrift

   | removes a lattice element and merges the two adjacent drift spaces
   
   |  **mergedrift** (SPLITPOS) removes an element located at SPLITPOS from the global lattice THERING
   |  surrounded by two DRIFT spaces. The resulting drift has Length L0 = L1 + LSPLIT + L2;
   |  Number of elements in THERING is thus reduced by 2
   
   | See also :func:`splitdrift`

.. py:function:: mvelem(elempos, dist)

   | Move an element
   
   | **mvelem(elempos, dist)** Move an element located at ELEMPOS in THERING
   |  surrounded by two DRIFT spaces
   
   |   0   < DIST  < LD move downstream
   |  -LU  < DIST  < 0  move upstream
   |   where LU and LD - lenths of
   |   upstream and downstrem drift drifts BEFORE!!! the move
   
   |  Number of elements in THERING and total length remain the same
   
   | See also :func:`splitdrift`, :func:`mergedrift`

.. py:function:: mvfield(dst,src,fieldnames)

   | Move fields from one structure to another
   
   | **[newdst,newsrc]=mvfield(dst,src,fieldnames)**
   |    Moves the selected fields from SRC to DST
   
   | DST:           Destination structure
   | SRC:           Source structure
   | FIELDNAMES:    Field names to be moved (Default: all fields from SRC)
   
   | NEWDST:        DST structure with fields added
   | NEWSRC:        SRC structure with fields removed

.. py:function:: rmelem0(lattice,elemindex)

   | removes elements of length 0 from the accelerator lattice
   |  **newlattice = rmelem0(lattice,elemindex)**
   |  **[newlattice, shiftedindex] = rmelem0(lattice,elemindex)**
   
   |  The number of elements in the modified lattice is
   |  reduced by length(ELEMINDEX)
   
   |  SHIFTEDINDEX points to elements in the NEWLATTICE that
   |  immediateley followed the removed elements in the original LATTICE
   
   | See also :func:`splitdrift`, :func:`mergedrift`

.. py:function:: setcellstruct(cellarray,...)

   | sets the field values of MATLAB cell array of structures
   
   |    Note that the calling syntax must be in the form of assignment:
   |    **cellarray = setcellstruct(cellarray,...)**
   |    MATLAB does not modify variables that only appear on the right
   |    hand side as arguments.
   
   |  Numeric data
   |  ---------------------------------------------------------
   |  **cellarray = setcellstruct(cellarray,'field',index,value,m,n)**
   | 	 Sets (M,N) element equal to VALUE when the field data is
   |    a matrix. The assigned VALUE may be
   |    1. Scalar numeric value to be written to all CELLARRAY elements
   |    2. Numeric array of the same length as INDEX array
   
   |  **cellarray = setcellstruct(cellarray,'field',index,value,m)** can be used
   |    for one dimensional vectors
   
   |  **cellarray = setcellstruct(cellarray,'field',index,value)** is the same as
   |    **setcellstruct(cellarray,'field',index,1,1)** if the field data
   |    is a scalar
   
   |  Character array
   |  --------------------------------------------------------------------
   |  CELLARRAY **setcellstruct(cellarray,'field',index,value)** is a MATLAB
   |    cell array of strings when specified fields contain strings.
   |    The assignment VALUE may be
   |    1. Character string,
   |    2. Character array with the number of rows matching the number of
   |        elements in INDEX array,
   |    3. Cell array of strings, with either one element or with the same
   |        length as index.
   
   | See also :func:`atsetfieldvalues`, :func:`getcellstruct`, :func:`findcells`

.. py:function:: setshift(elemindex, dx, dy)

   | sets the misalignment vectors T1, T2 for elements
   
   |  **setshift(elemindex, dx, dy)** sets the entrance and exit misalignment vectors
   |   of one element or a group of elements in the globally defined lattice THERING.
   
   |   DX, DY are displacements of the ELEMENT
   |   so the coordinate transformation on the particle at entrance is
   | 	X  ->  X-DX
   |    Y  ->  Y-DY
   |   The elements to be modified are given by ELEMINDEX
   | 	Previous stored values are overwritten.
   
   | See also :func:`setshift`

.. py:function:: settags(lattice, index, tag)

   | sets the 'Tag' field in AT lattice elements
   |  **lattice = settags(lattice, index, tag)**
   |    INDEX can be integer AT index or a string famly name
   |    TAG is a string tag or a cell array of strings
   |  **lattice = settags(lattice, index, tag, 'append')**
   |    appends to existing tags

.. py:function:: settilt(elemindex, psi)

   | sets the entrance and exit misalignment matrixes
   |  of an element or a group of elements in THERING
   |  Previously stored values are overwritten.
   
   |  **settilt(elemindex, psi)**
   |  ELEMINDEX contains indexes of elements to be rotated
   |  PSI - angle(s) of rotation in RADIANS
   |    POSITIVE PSI corresponds to a CORKSCREW (right)
   |    rotation of the ELEMENT.
   |    (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
   |    The misalgnment matrixes are stored in fields R1 and R2
   |    R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
   |    R2 = R1'
   
   |   NOTES
   |   1. This function is deprecated. Use atsettilt instead
   
   | See also :func:`setshift`, :func:`mksrollmat`

.. py:function:: splitdrift(driftpos, split) inserts a marker (zero-length)

   | inserts an element into a drift space
   
   |  **splitdrift(driftpos, split) inserts a marker (zero-length)** element
   |    at distance SPLIT ( 0 < SPLIT < 1) into a drift space
   |    located at DRIFTPOS in THERING
   
   |  **splitdrift(driftpos, split, elemstruccture) inserts a marker (zero-length)** element
   |    at distance SPLIT ( 0 < SPLIT < 1) into a drift space
   |    located at DRIFTPOS in THERING
   
   |  Number of elements in the RING is thus increased by 2
   |  SPLIT (controls the position of the split
   |  L1 = L0*SPLIT
   |  L2 = L0(1-SPLIT)
   |   where L0 is the length of the original DRIFT
   
   | See also :func:`mergedrift`

.. py:function:: symplectify

   | symplectify makes a matrix more symplectic
   | follow Healy algorithm as described by McKay
   | BNL-75461-2006-CP

