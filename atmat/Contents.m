% Accelerator toolbox
% Version 2.1.002 (atcollab) 01-Aug-2021
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
%*Element creation*
%       atbaselem                - Create an AT element structure + various checks
%       atelem                   - Makes a new AT element structure from another element,
%       atringparam              - Creates a RingParameter Element which should go at the beginning of the ring
%       atdrift                  - Creates a drift space element with Class 'Drift'
%       atquadrupole             - Creates a quadrupole element with Class 'Quadrupole'
%       atrbend                  - Creates a rectangular bending magnet element with class 'Bend'
%       atsbend                  - Creates a sector bending magnet element with class 'Bend'
%       atsextupole              - Creates a sextupole element with class 'Sextupole'
%       atskewquad               - Creates a skew quadrupole element with Class 'Multipole'
%       atmultipole              - Creates a multipole element
%       atrfcavity               - Creates an rfcavity element with Class 'RFCavity'
%       atmarker                 - Creates a marker space element
%       atmonitor                - Creates a Beam Position Monitor element with Class 'Monitor'
%       ataperture               - Creates a aperture element
%       atcorrector              - Creates a drift space element with class 'Corrector'
%       atidtable                - Creates an ID element
%       atwiggler                - Creates a wiggler
%       atdampMatElem            - Creates an element that applies the global damping matrix
%       atsolenoid               - Creates a new solenoid element with Class 'Solenoid'
%       atthinmultipole          - Creates a thin multipole element
%       atM66                    - Create an element applying an arbitrary 6x6 transfer matrix
%       atQuantDiff              - Creates a quantum diffusion element
%
%*Element manipulation*
%       isatelem                 - Tests if an input argument is a valid AT element
%       atguessclass             - Tries to determine the class of an element
%       atshiftelem              - Set new displacement parameters
%       attiltelem               - Sets new rotation parameters
%
%*Lattice manipulation*
%       isatlattice              - Tests if an input argument is a valid AT lattice
%   Access elements
%       atindex                  - Extracts the information about element families and
%       atgetcells               - Performs a search on MATLAB cell arrays of structures
%       atgetfieldvalues         - Retrieves the field values AT cell array of elements
%       atsetfieldvalues         - Sets the field values of MATLAB cell array of structures
%   Insert elements
%       atinsertelems            - Insert elements at given locations in a line
%       atdivelem                - Divide an element into pieces
%       atsplitelem              - Creates a line by inserting one or more elements into a base element
%       insertindrift            - Inserts one or more elements into a drift element
%       atsbreak                 - Insert markers at given s positions in a lattice
%   Join elements
%       atreduce                 - Remove useless elements from an AT structure
%       combinebypassmethod      - Combines adjacent elements that have the same specified pass method
%       combinelinear45          - Combines adjacent  elements that use 4-by-5 PassMethods
%   Other
%       atloadfielderrs          - Will load a field error structure into a ring
%       atsetRFCavity            - Set the RF Cavity with the passmethod RFCavityPass
%       atsetshift               - Sets the misalignment vectors
%       atsettilt                - Sets the entrance and exit rotation matrices
%       settags                  - Sets the 'Tag' field in AT lattice elements
%       findtags                 - Looks for string matches in 'Tag' field of AT lattice elements
%       mvelem                   - Move an element
%       mvfield                  - Move fields from one structure to another
%
%*Loading and Saving lattices*
%   Binary files
%
% Lattices can be saved as binary mat-files using the standard load and save commands
%
%   Text Files
%       atwritem                 - Creates a .m file to store an AT structure
%       atwritepy                - Creates pyAT lattice from a Matlab lattice
%
%*Linear optics*
%   Closed orbit
%       findorbit                - Find the closed orbit
%       findorbit4               - Finds closed orbit in the 4-d transverse phase
%       findorbit6               - Finds closed orbit in the full 6-d phase space
%       findsyncorbit            - Finds closed orbit, synchronous with the RF cavity
%   Transfer matrices
%       findm44                  - Numerically finds the 4x4 transfer matrix of an accelerator lattice
%       findm66                  - Numerically finds the 6x6 transfer matrix of an accelerator lattice
%       findelemm44              - Numerically finds the 4x4 transfer matrix of an element
%       findelemm66              - Numerically finds the 6x6 transfer matrix of an element
%   Optical functions
%       atlinopt2                - Performs the linear analysis of UNCOUPLED lattices
%       atlinopt4                - Performs the 4D linear analysis of COUPLED lattices
%       atlinopt6                - Performs linear analysis of the lattice
%       beam22                   - Computes the beam matrix from the 1-turn transfer matrix
%       beam44                   - Computes the coupled beam matrices
%
%*Radiation*
%       check_radiation          - Check the radiation state of a ring
%       atenergy                 - Gets the lattice energy
%       atgetU0                  - Computes Energy loss per turn in eV 
%       atdampingrates           - Find tunes and damping rates from one map matrix with radiation
%       atradon                  - Switches RF and radiation on
%       atradoff                 - Switches radiation and cavity off
%       quantumDiff              - Compute the radiation-diffusion matrix
%       ohmienvelope             - Calculates equilibrium beam envelope in a
%       DipoleRadiation          - Compute the radiation integrals in dipoles
%       WigglerRadiation         - Compute the radiation integrals in wigglers
%
%*Parameter summary*
%       atx                      - Computes and displays global information
%       atsummary                - Print out the parameters of the current AT lattice
%       ringpara                 - Calculates various ring parameters
%
%*Physics*
%       symplectify              - Makes a matrix more symplectic
