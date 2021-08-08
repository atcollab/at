% Accelerator toolbox
% Version 2.0 (atcollab) 01-Aug-2021
%
% The Accelerator Toolbox was originally created by Andrei Terebilo.
% Development is now continued by a multi-laboratory collaboration, atcollab
%
% The integrators used for tracking particles may be compiled with
% the command atmexall.
%
% For getting started, one may look at the examples in atmat/atdemos.
% Example lattice files are located in machine_data.
%
% *Element creation*
%%
%
% <matlab:help('atbaselem') atbaselem>                      - Create an AT element structure + various checks
% <matlab:help('atelem') atelem>                         - Makes a new AT element structure from another element,
% <matlab:help('atringparam') atringparam>                    - Creates a RingParameter Element which should go at the beginning of the ring
% <matlab:help('atdrift') atdrift>                        - Creates a drift space element with Class 'Drift'
% <matlab:help('atquadrupole') atquadrupole>                   - Creates a quadrupole element with Class 'Quadrupole'
% <matlab:help('atrbend') atrbend>                        - Creates a rectangular bending magnet element with class 'Bend'
% <matlab:help('atsbend') atsbend>                        - Creates a sector bending magnet element with class 'Bend'
% <matlab:help('atsextupole') atsextupole>                    - Creates a sextupole element with class 'Sextupole'
% <matlab:help('atskewquad') atskewquad>                     - Creates a skew quadrupole element with Class 'Multipole'
% <matlab:help('atmultipole') atmultipole>                    - Creates a multipole element
% <matlab:help('atrfcavity') atrfcavity>                     - Creates an rfcavity element with Class 'RFCavity'
% <matlab:help('atmarker') atmarker>                       - Creates a marker space element
% <matlab:help('atmonitor') atmonitor>                      - Creates a Beam Position Monitor element with Class 'Monitor'
% <matlab:help('ataperture') ataperture>                     - Creates a aperture element
% <matlab:help('atcorrector') atcorrector>                    - Creates a drift space element with class 'Corrector'
% <matlab:help('atidtable') atidtable>                      - (FAMNAME,Nslice,filename,Energy,method)
% <matlab:help('atwiggler') atwiggler>                      - Creates a wiggler
% <matlab:help('atdampMatElem') atdampMatElem>                  - Creates an element that applies the global damping matrix
% <matlab:help('atsolenoid') atsolenoid>                     - Creates a new solenoid element with Class 'Solenoid'
% <matlab:help('atthinmultipole') atthinmultipole>                - Creates a thin multipole element
% <matlab:help('atM66') atM66>                          - (FAMNAME,M66,PASSMETHOD)
% <matlab:help('atQuantDiff') atQuantDiff>                    - Creates a quantum diffusion element
%
% *Element manipulation*
%
% <matlab:help('isatelem') isatelem>                       - Tests if an input argument is a valid AT element
% <matlab:help('atguessclass') atguessclass>                   - Tries to determine the class of an element
% <matlab:help('atshiftelem') atshiftelem>                    - Set new displacement parameters
% <matlab:help('attiltelem') attiltelem>                     - Sets new rotation parameters
%
% *Lattice manipulation*
%
% <matlab:help('isatlattice') isatlattice>                    - Tests if an input argument is a valid AT lattice
%      Access elements
% <matlab:help('atindex') atindex>                        - Extracts the information about element families and
% <matlab:help('atgetcells') atgetcells>                     - Performs a search on MATLAB cell arrays of structures
% <matlab:help('atgetfieldvalues') atgetfieldvalues>               - Retrieves the field values AT cell array of elements
% <matlab:help('atsetfieldvalues') atsetfieldvalues>               - Sets the field values of MATLAB cell array of structures
%      Insert elements
% <matlab:help('atinsertelems') atinsertelems>                  - Insert elements at given locations in a line
% <matlab:help('atdivelem') atdivelem>                      - LINE=ATDIVELEM(ELEM,FRAC) divide an element into pieces
% <matlab:help('atsplitelem') atsplitelem>                    - Creates a line by inserting one or more elements into a base element
% <matlab:help('insertindrift') insertindrift>                  - Inserts one or more elements into a drift element
% <matlab:help('atsbreak') atsbreak>                       - Insert markers at given s positions in a lattice
%      Join elements
% <matlab:help('atreduce') atreduce>                       - Remove useless elements from an AT structure
% <matlab:help('mergedrift') mergedrift>                     - Removes a lattice element and merges the two adjacent drift spaces
% <matlab:help('combinebypassmethod') combinebypassmethod>            - Combines adjacent elements that have the same specified pass method
% <matlab:help('combinelinear45') combinelinear45>                - Combines adjacent  elements that use 4-by-5 PassMethods
%      Other
% <matlab:help('atloadfielderrs') atloadfielderrs>                - Will load a field error structure into a ring
% <matlab:help('atsetRFCavity') atsetRFCavity>                  - - sets the RF Cavity with the passmethod RFCavityPass
% <matlab:help('atsetshift') atsetshift>                     - Sets the misalignment vectors
% <matlab:help('atsettilt') atsettilt>                      - Sets the entrance and exit rotation matrices
% <matlab:help('settags') settags>                        - Sets the 'Tag' field in AT lattice elements
% <matlab:help('findtags') findtags>                       - Looks for string matches in 'Tag' field of AT lattice elements
% <matlab:help('mvelem') mvelem>                         - (ELEMPOS, DIST) moves an element  located at ELEMPOS in THERING
% <matlab:help('mvfield') mvfield>                        - Move fields from one structure to another
% <matlab:help('reverse') reverse>                        - Reverses the order of elements in a one-dimensional MATLAB ARRAY
% <matlab:help('splitdrift') splitdrift>                     - Inserts an element into a drift space
%
% *Loading and Saving lattices*
%
% <matlab:help('atwritem') atwritem>                       - Creates a .m file to store an AT structure
% <matlab:help('atwritepy') atwritepy>                      - Creates a .m file to store an AT structure
%
% *Linear optics*
%
% <matlab:help('atlinopt2') atlinopt2>                      - Performs the linear analysis of UNCOUPLED lattices
% <matlab:help('atlinopt4') atlinopt4>                      - Performs the 4D linear analysis of COUPLED lattices
% <matlab:help('atlinopt6') atlinopt6>                      - Performs linear analysis of the lattice
% <matlab:help('beam22') beam22>                         - Computes the beam matrix from the 1-turn transfer matrix
% <matlab:help('beam44') beam44>                         - Computes the coupled beam matrices
% <matlab:help('findm44') findm44>                        - Numerically finds the 4x4 transfer matrix of an accelerator lattice
% <matlab:help('findm66') findm66>                        - Numerically finds the 6x6 transfer matrix of an accelerator lattice
% <matlab:help('findelemm44') findelemm44>                    - Numerically finds the 4x4 transfer matrix of an element
% <matlab:help('findelemm66') findelemm66>                    - Numerically finds the 6x6 transfer matrix of an element
%
% *Physics*
%
% <matlab:help('symplectify') symplectify>                    - Makes a matrix more symplectic
