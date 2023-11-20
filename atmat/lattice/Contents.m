% LATTICE
% See also 
% 
%   Contents file for LATTICE and its subfolders.
%   
%   LATTICE
%   at2py                      - ELSTR=AT2PY(ELEM) convert AT element tp pyat
%   at2str                     - Makes the string representation of an AT element
%   ataddmpolecomppoly         - Adds a multipole component to an existing polynomial,
%   ataddmpoleerrors           - Ataddrandmpole adds a random multipole component to all elements of type
%   atCheckRingProperties      - Get the ring properties if existing
%   atdivelem                  - Divide an element into pieces
%   atelem                     - Makes a new AT element structure from another element,
%   atfastring                 - Generate simplified AT structures
%   atfitchrom                 - Fit chromaticites by scaling 2 sextupole families
%   atfittune                  - Fit linear tunes by scaling 2 quadrupole families
%   atgetcells                 - Performs a search on MATLAB cell arrays of structures
%   atgetfieldvalues           - Retrieves the field values AT cell array of elements
%   atGetRingProperties        - Get the ring properties
%   atguessclass               - Tries to determine the class of an element
%   atindex                    - Extracts the information about element families and
%   atinsertelems              - Insert elements at given locations in a line
%   atloadfielderrs            - Will load a field error structure into a ring
%   atloadlattice              - Load a lattice from a file
%   atlocateparam              - Private function. Locate the RingParam element
%   atmaincavities             - Get the fundamental mode cavities
%   atmakefielderrstruct       - MAKERNDFIELDERRS will create a field error data structure
%   atparamscan                - Private function. Updates the RingParam element
%   atparticle                 - Particle definition for AT
%   atreduce                   - Remove useless elements from an AT structure
%   atrotatelattice            - Circularly shift the lattice elements
%   atsbreak                   - Insert markers at given s positions in a lattice
%   atsetfieldvalues           - Sets the field values of MATLAB cell array of structures
%   atsetRFCavity              - Set the RF Cavity with the passmethod RFCavityPass
%   atSetRingProperties        - Add or modify properties of the lattice
%   atsetshift                 - Sets the misalignment vectors
%   atsettilt                  - Sets the entrance and exit rotation matrices
%   atshiftelem                - Set new displacement parameters
%   atsimplering               - Creates a "simple ring"
%   atsplitelem                - Creates a line by inserting one or more elements into a base element
%   attiltelem                 - Sets new rotation parameters
%   atwritem                   - Creates a .m file to store an AT structure
%   atwritepy                  - Creates pyAT lattice from a Matlab lattice
%   buildlat                   - Places elements from FAMLIST into cell array THERING
%   combinebypassmethod        - Combines adjacent elements that have the same specified pass method
%   combinelinear45            - Combines adjacent  elements that use 4-by-5 PassMethods
%   entrancefields             - () Return the list of field names affecting the element entrance
%   exitfields                 - () Return the list of field names affecting the element exit
%   findcells                  - Performs a search on MATLAB cell arrays of structures
%   findtags                   - Looks for string matches in 'Tag' field of AT lattice elements
%   getcellstruct              - Retrieves the field values MATLAB cell array of structures
%   insertelem0                - - quick and dirty:
%   insertindrift              - Inserts one or more elements into a drift element
%   isatelem                   - Tests if an input argument is a valid AT element
%   isatlattice                - Tests if an input argument is a valid AT lattice
%   mergedrift                 - Removes a lattice element and merges the two adjacent drift spaces
%   mvelem                     - Move an element
%   mvfield                    - Move fields from one structure to another
%   rmelem0                    - Removes elements of length 0 from the accelerator lattice
%   setcellstruct              - Sets the field values of MATLAB cell array of structures
%   setshift                   - Sets the misalignment vectors T1, T2 for elements
%   settags                    - Sets the 'Tag' field in AT lattice elements
%   settilt                    - Sets the entrance and exit misalignment matrixes
%   splitdrift                 - Inserts an element into a drift space
%   symplectify                - Makes a matrix more symplectic
%   
%   LATTICE/CONVERTERS
%   readmad                    - Reads the file output of MAD commands
%   
%   LATTICE/CONVERTERS/AT2ELEGANT
%   AT_2_Elegant               - This functions converts the AT lattice AT_ring in elegant form
%   
%   LATTICE/CONVERTERS/AT2G4BL
%   ATtoG4BL                   - Function [outtext]=ATtoG4BL(P_0,particle,folder)
%   
%   LATTICE/CONVERTERS/AT2MAD8
%   AT_2_mad8                  - Function [elelat,def,lines]=AT_2_mad8(AT_ring,linename)
%   
%   LATTICE/CONVERTERS/AT2MADX
%   AT_2_madX                  - Function [elelat,defs,lines]=AT_2_madX(AT_ring,linename)
%   
%   LATTICE/CONVERTERS/AT2OPA
%   AT_2_OPA                   - Function AT_2_OPA(AT_ring,linename)
%   
%   LATTICE/CONVERTERS/ELEGANT2AT
%   ele2at_run_me              - Test_elegant_converter
%   elegant2at                 - Function elegant2at(elegantlattice,E0,outfilename)
%   ParseAtributesELEGANT_2_AT - Determines atribute and sets field in sxs{i} structure AT
%   
%   LATTICE/CONVERTERS/MAD82MADX
%   mad8TOmadx                 - Converts mad8 sequence files to madX
%   
%   LATTICE/CONVERTERS/MADX2AT
%   atfrommadx                 - Function atfrommadx(seqfilemadX,E0,outfilename)
%   buildATLattice             - Given a list (cell array) of elements with specified field Spos (center of element (madx default)) in a
%   ParseAtributesMADX_2_AT    - Determines atribute and sets field in sxs{i} structure AT
%   reshapeToCellArray         - If CEL_CEL is a cell array of structures and cell arrays it converts it a
%   
%   LATTICE/CONVERTERS/MADX2AT/EXAMPLES
%   convertMADXtoATExample     - Simple lattice test (uncomment to run)
%   
%   LATTICE/CONVERTERS/MADX2G4BL
%   madx2g4bl                  - Function [outtext]=madx2g4bl(P_0,particle,folder)
%   
%   LATTICE/ELEMENT_CREATION
%   ataperture                 - Creates a aperture element
%   atbaselem                  - Create an AT element structure + various checks
%   atcorrector                - Creates a drift space element with class 'Corrector'
%   atdampMatElem              - Creates an element that applies the global damping matrix
%   atdrift                    - Creates a drift space element with Class 'Drift'
%   atenergyloss               - Creates an energy loss element
%   atidtable                  - Creates an ID element
%   atinsertiondevicekickmap   - Creates an insertion device kick-map element
%   atM66                      - Create an element applying an arbitrary 6x6 transfer matrix
%   atM66Tijk                  - ATM66(FAMNAME,M66,Tijk,PASSMETHOD)
%   atmarker                   - Creates a marker space element
%   atmonitor                  - Creates a Beam Position Monitor element with Class 'Monitor'
%   atmultipole                - Creates a multipole element
%   atquadrupole               - Creates a quadrupole element with Class 'Quadrupole'
%   atQuantDiff                - Creates a quantum diffusion element
%   atrbend                    - Creates a rectangular bending magnet element with class 'Bend'
%   atrbendtune                - Set X0ref and RefDZ for rectangular bending magnets
%   atrfcavity                 - Creates an rfcavity element with Class 'RFCavity'
%   atringparam                - Creates a RingParameter Element which should go at the beginning of the ring
%   atsbend                    - Creates a sector bending magnet element with class 'Bend'
%   atsextupole                - Creates a sextupole element with class 'Sextupole'
%   atSimpleQuantDiff          - SimpleQuantDiff creates a simple quantum difusion element
%   atskewquad                 - Creates a skew quadrupole element with Class 'Multipole'
%   atsolenoid                 - Creates a new solenoid element with Class 'Solenoid'
%   atthinmultipole            - Creates a thin multipole element
%   atvariablemultipole        - Creates a variable thin multipole element
%   atwiggler                  - Creates a wiggler
%   
%   LATTICE/ELEMENT_CREATION/PRIVATE
%   decodeatargs               - Separates arguments and resources
%   
%   LATTICE/PARAMGROUP
%   atparamgroup               - PARAMETER GROUP in AT is a general way
%   mkparamgroup               - Simplifies creation of AT parameter groups
%   restoreparamgroup          - Restores the values of multiple physical
%   saveparamgroup             - Saves the values of multiple physical
%   setparamgroup              - Modifies a group of parameters
%   
%   LATTICE/SURVEY
%   atgeometry                 - Computes the 2-D position of all elements (no vertical bend)
%   atgeometry3                - Computes the 3-D position of all elements
%    
%   This file was generated by updateContents.m on 06 Nov 2023 at 14:58:12.
