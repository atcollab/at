% Accelerator Toolbox
% Version 2.2-dev.354 (atcollab) 01-Feb-2022
%
% The Accelerator Toolbox was originally created by Andrei Terebilo.
% Development is now continued by a multi-laboratory collaboration, <a href="matlab:web('https://github.com/atcollab')">atcollab</a>
%
% <a href="matlab:web('https://atcollab.github.io/at/')">AT web site</a>
%
% The integrators used for tracking particles may be compiled with
% the command atmexall.
%
% For getting started, one may look at the examples in atmat/atdemos.
% Example lattice files are located in machine_data.
%
%*Element creation*
%
%       <a href="matlab:help atbaselem">atbaselem</a>                - Create an AT element structure + various checks
%       <a href="matlab:help atelem">atelem</a>                   - Makes a new AT element structure from another element,
%       <a href="matlab:help atringparam">atringparam</a>              - Creates a RingParameter Element which should go at the beginning of the ring
%       <a href="matlab:help atdrift">atdrift</a>                  - Creates a drift space element with Class 'Drift'
%       <a href="matlab:help atquadrupole">atquadrupole</a>             - Creates a quadrupole element with Class 'Quadrupole'
%       <a href="matlab:help atrbend">atrbend</a>                  - Creates a rectangular bending magnet element with class 'Bend'
%       <a href="matlab:help atsbend">atsbend</a>                  - Creates a sector bending magnet element with class 'Bend'
%       <a href="matlab:help atsextupole">atsextupole</a>              - Creates a sextupole element with class 'Sextupole'
%       <a href="matlab:help atskewquad">atskewquad</a>               - Creates a skew quadrupole element with Class 'Multipole'
%       <a href="matlab:help atmultipole">atmultipole</a>              - Creates a multipole element
%       <a href="matlab:help atrfcavity">atrfcavity</a>               - Creates an rfcavity element with Class 'RFCavity'
%       <a href="matlab:help atmarker">atmarker</a>                 - Creates a marker space element
%       <a href="matlab:help atmonitor">atmonitor</a>                - Creates a Beam Position Monitor element with Class 'Monitor'
%       <a href="matlab:help ataperture">ataperture</a>               - Creates a aperture element
%       <a href="matlab:help atcorrector">atcorrector</a>              - Creates a drift space element with class 'Corrector'
%       <a href="matlab:help atidtable">atidtable</a>                - Creates an ID element
%       <a href="matlab:help atwiggler">atwiggler</a>                - Creates a wiggler
%       <a href="matlab:help atdampMatElem">atdampMatElem</a>            - Creates an element that applies the global damping matrix
%       <a href="matlab:help atsolenoid">atsolenoid</a>               - Creates a new solenoid element with Class 'Solenoid'
%       <a href="matlab:help atthinmultipole">atthinmultipole</a>          - Creates a thin multipole element
%       <a href="matlab:help atM66">atM66</a>                    - Create an element applying an arbitrary 6x6 transfer matrix
%       <a href="matlab:help atQuantDiff">atQuantDiff</a>              - Creates a quantum diffusion element
%
%*Element manipulation*
%
%       <a href="matlab:help isatelem">isatelem</a>                 - Tests if an input argument is a valid AT element
%       <a href="matlab:help atguessclass">atguessclass</a>             - Tries to determine the class of an element
%       <a href="matlab:help atshiftelem">atshiftelem</a>              - Set new displacement parameters
%       <a href="matlab:help attiltelem">attiltelem</a>               - Sets new rotation parameters
%
%*Lattice manipulation*
%
%       <a href="matlab:help isatlattice">isatlattice</a>              - Tests if an input argument is a valid AT lattice
%
%   Global lattice parameters
%       The global lattice properties 'FamName', 'Energy', 'Periodicity', 'Particle', 'HarmNumber' are stored in the RingParam lattice element.
%       The following functions gives an easy access to them:
%
%       <a href="matlab:help atGetRingProperties">atGetRingProperties</a>      - Get the ring properties
%       <a href="matlab:help atSetRingProperties">atSetRingProperties</a>      - Add or modify properties of the lattice
%       <a href="matlab:help atenergy">atenergy</a>                 - Gets the lattice energy
%
%   Access elements
%       <a href="matlab:help atindex">atindex</a>                  - Extracts the information about element families and
%       <a href="matlab:help atgetcells">atgetcells</a>               - Performs a search on MATLAB cell arrays of structures
%       <a href="matlab:help atgetfieldvalues">atgetfieldvalues</a>         - Retrieves the field values AT cell array of elements
%       <a href="matlab:help atsetfieldvalues">atsetfieldvalues</a>         - Sets the field values of MATLAB cell array of structures
%
%   Insert elements
%       <a href="matlab:help atinsertelems">atinsertelems</a>            - Insert elements at given locations in a line
%       <a href="matlab:help atdivelem">atdivelem</a>                - Divide an element into pieces
%       <a href="matlab:help atsplitelem">atsplitelem</a>              - Creates a line by inserting one or more elements into a base element
%       <a href="matlab:help insertindrift">insertindrift</a>            - Inserts one or more elements into a drift element
%       <a href="matlab:help atsbreak">atsbreak</a>                 - Insert markers at given s positions in a lattice
%
%   Join elements
%       <a href="matlab:help atreduce">atreduce</a>                 - Remove useless elements from an AT structure
%       <a href="matlab:help combinebypassmethod">combinebypassmethod</a>      - Combines adjacent elements that have the same specified pass method
%       <a href="matlab:help combinelinear45">combinelinear45</a>          - Combines adjacent  elements that use 4-by-5 PassMethods
%
%   Other
%       <a href="matlab:help atloadfielderrs">atloadfielderrs</a>          - Will load a field error structure into a ring
%       <a href="matlab:help atsetRFCavity">atsetRFCavity</a>            - Set the RF Cavity with the passmethod RFCavityPass
%       <a href="matlab:help atsetshift">atsetshift</a>               - Sets the misalignment vectors
%       <a href="matlab:help atsettilt">atsettilt</a>                - Sets the entrance and exit rotation matrices
%       <a href="matlab:help settags">settags</a>                  - Sets the 'Tag' field in AT lattice elements
%       <a href="matlab:help findtags">findtags</a>                 - Looks for string matches in 'Tag' field of AT lattice elements
%       <a href="matlab:help mvelem">mvelem</a>                   - Move an element
%       <a href="matlab:help mvfield">mvfield</a>                  - Move fields from one structure to another
%
%*Loading and Saving lattices*
%
%   Binary files
%       Lattices can be saved as binary mat-files using the standard load and save commands
%
%   Text Files
%       <a href="matlab:help atwritem">atwritem</a>                 - Creates a .m file to store an AT structure
%       <a href="matlab:help atwritepy">atwritepy</a>                - Creates pyAT lattice from a Matlab lattice
%
%*Linear optics*
%
%   Closed orbit
%       <a href="matlab:help findorbit">findorbit</a>                - Find the closed orbit
%       <a href="matlab:help findorbit4">findorbit4</a>               - Finds closed orbit in the 4-d transverse phase
%       <a href="matlab:help findorbit6">findorbit6</a>               - Finds closed orbit in the full 6-d phase space
%       <a href="matlab:help findsyncorbit">findsyncorbit</a>            - Finds closed orbit, synchronous with the RF cavity
%
%   Transfer matrices
%       <a href="matlab:help findm44">findm44</a>                  - Numerically finds the 4x4 transfer matrix of an accelerator lattice
%       <a href="matlab:help findm66">findm66</a>                  - Numerically finds the 6x6 transfer matrix of an accelerator lattice
%       <a href="matlab:help findelemm44">findelemm44</a>              - Numerically finds the 4x4 transfer matrix of an element
%       <a href="matlab:help findelemm66">findelemm66</a>              - Numerically finds the 6x6 transfer matrix of an element
%
%   Optical functions
%       <a href="matlab:help atlinopt2">atlinopt2</a>                - Performs the linear analysis of UNCOUPLED lattices
%       <a href="matlab:help atlinopt4">atlinopt4</a>                - Performs the 4D linear analysis of COUPLED lattices
%       <a href="matlab:help atlinopt6">atlinopt6</a>                - Performs linear analysis of the lattice
%       <a href="matlab:help beam22">beam22</a>                   - Computes the beam matrix from the 1-turn transfer matrix
%       <a href="matlab:help beam44">beam44</a>                   - Computes the coupled beam matrices
%
%*Radiation*
%
%       <a href="matlab:help check_radiation">check_radiation</a>          - Check the radiation state of a ring
%       <a href="matlab:help atenergy">atenergy</a>                 - Gets the lattice energy
%       <a href="matlab:help atgetU0">atgetU0</a>                  - Computes Energy loss per turn in eV 
%       <a href="matlab:help atdampingrates">atdampingrates</a>           - Find tunes and damping rates from one map matrix with radiation
%       <a href="matlab:help atradon">atradon</a>                  - Switches RF and radiation on
%       <a href="matlab:help atradoff">atradoff</a>                 - Switches radiation and cavity off
%       <a href="matlab:help quantumDiff">quantumDiff</a>              - Compute the radiation-diffusion matrix
%       <a href="matlab:help ohmienvelope">ohmienvelope</a>             - Calculates equilibrium beam envelope in a
%       <a href="matlab:help DipoleRadiation">DipoleRadiation</a>          - Compute the radiation integrals in dipoles
%       <a href="matlab:help WigglerRadiation">WigglerRadiation</a>         - Compute the radiation integrals in wigglers
%
%*Parameter summary*
%
%       <a href="matlab:help atx">atx</a>                      - Computes and displays global information
%       <a href="matlab:help atsummary">atsummary</a>                - Print out the parameters of the current AT lattice
%       <a href="matlab:help ringpara">ringpara</a>                 - Calculates various ring parameters
%
%*Physics*
%
%       <a href="matlab:help symplectify">symplectify</a>              - Makes a matrix more symplectic
%
%<a href="matlab:web('/Applications/MATLAB_R2021b.app/help/3ptoolbox/atacceleratortoolbox/doc/AT_page.html')">See documentation for AT</a>
