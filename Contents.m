% Accelerator Toolbox
% Version 2.5.0 (atcollab) 05-Nov-2023
% 
% AT
% See also 
% 
%   Contents file for AT and its subfolders.
%   
%   AT/ATINTEGRATORS
%   aperture                            - Creates a aperture element in old a AT version (Obsolete)
%   BndStrMPoleSymplectic4Pass          - .m Help file for BndStrMPoleSymplectic4Pass.c
%   corrector                           - Creates a corrector element in old a AT version (Obsolete)
%   drift                               - Creates a drift element in old a AT version (Obsolete)
%   DriftPass                           - - Integrator for Drift spaces
%   ExactHamiltonianPass                - .m Help file for ExactHamiltonianPass.c
%   GWigSymplecticPass                  - #ok<STOUT>
%   GWigSymplecticRadPass               - #ok<STOUT>
%   hmonitor                            - Creates a horizontal monitor element in old AT version (Obsolete)
%   idtable                             - (fname, Nslice, filename, Energy, method)
%   idtable_global                      - IDTABLE Creates a RADIA-Map based element
%   marker                              - Creates a marker element in old AT version (Obsolete)
%   mChangeRefPass                      - Change the reference energy by scaling
%   mDriftPass                          - Example of pass method in matlab
%   mhdrload_bis                        - Loads data from an ASCII file containing multiple text
%   mIdentityPass                       - MDRIFTPASS - example of pass method in matlab
%   monitor                             - Creates a monitor element in old AT versions (Obsolete)
%   multipole                           - Creates a thick multipole element in old AT versions (Obsolete)
%   passmethodlist                      - Utility function for MATLAB Compiler
%   passmethods                         - Returns a list of available AT passmethod functions in
%   quadrupole                          - Creates a quadrupole element in old AT version (Obsolete)
%   rbend                               - Creates a rectangular bend in old AT version (Obsolete)
%   rbend2                              - Creates rectangular bend with dipole fringe field  in old AT versions (Obsolete)
%   rbend3                              - Creates a rectangular bend with different fringe fields at entrance
%   rfcavity                            - Creates a RF cavity in older AT version
%   sbend                               - Creates a sector bend element in old AT versions (Obsolete)
%   sextupole                           - Creates a sextupole element in old AT versions (Obsolete)
%   skewquad                            - Creates a skewquad element (alias to multipole) in old AT version (Obsolete)
%   solenoid                            - Creates a solenoid element in old AT versions (Obsolete)
%   vmonitor                            - Creates a vertical monitor element in old AT version (Obsolete)
%   wiggler                             - (fname, Ltot, Lw, Bmax, Nstep, Nmeth, By, Bx, method,energy)
%   
%   AT/ATMAT
%   at                                  - Accelerator Toolbox
%   atcavityoff                         - Switches RF cavities off
%   atcavityon                          - ATRADON switches RF cavities on
%   atdiag                              - Tests AT intallation
%   atdisplay                           - Checks the verbosity level in the global variable GLOBVAL
%   athelp                              - Generate the list of Accelerator Toolbox functions
%   atm2html                            - MAKEDOC_HTML - Generate new MML, SOLEIL and AT HTML help files
%   atmexall                            - Build all AT platform dependent mex-files from C-sources
%   atpath                              - Adds the AT directories to the MATLAB search path
%   atroot                              - Returns Accelerator Toolbox root directory
%   atupdateContents                    - Updates Contents.m of AT directories
%   getContents                         - Get the contents of a specified directory
%   isOctave                            - Check if running Octave
%   updateContents                      - Create a Contents.m file including subdirectories
%   
%   AT/ATMAT/ATDEMOS
%   
%   AT/ATMAT/ATDEMOS/ATFASTRING
%   testatfastring                      - Demo of using atfastring for fast tracking
%   
%   AT/ATMAT/ATDEMOS/ATMATCHEXAMPLES/BUMP
%   run_BumpFit                         - Fit a bump using correctors
%   
%   AT/ATMAT/ATDEMOS/ATMATCHEXAMPLES/EXAMPLEATMATCH
%   betx                                - Get value of betx for  Seq(indx)
%   bety                                - Get value of bety for  Seq(indx)
%   dispx                               - Get value of horizontal dispersion for  Seq(indx)
%   getDispersion                       - Function [dx,dy]=getDispersion(THERING,bpmindx)
%   mux
%   runtotest_atmatch                   - Macro match dba test lattice beta functions and dispersion using
%   VaryQuadFam
%   
%   AT/ATMAT/ATDEMOS/ATMATCHEXAMPLES/MATCHCHROMATICITY
%   run_matchChromaticity               - Match chrmaticity
%   
%   AT/ATMAT/ATDEMOS/GUI
%   atslider                            - Is an example of a GUI control of multiple parameters in THERING
%   demoknob                            - Illustrates the use of MATLAB GUI controls with AT
%   
%   AT/ATMAT/ATDEMOS/IDMODELING
%   add_ID                              - Adds IDelem at the beginning and subtracts the length from
%   trackWithU35KickMap                 - Read in ESRF lattice
%   
%   AT/ATMAT/ATDEMOS/OPTICSANDBEAMSIZES
%   latticedemo                         - Self-running tutorial
%   linoptdemo                          - Script illustrates the use of LINOPT
%   ohmienvelopedemo                    - Illustrates the use of OHMIENVELOPE function
%   
%   AT/ATMAT/ATDEMOS/RESPONSEMATRIX
%   findrespmdemo                       - Response matrix demo
%   
%   AT/ATMAT/ATDEMOS/TRACKING
%   generateTrackData                   - Give some initial coordinates.  Track through sample lattices
%   test_openmp                         - Test the perfoemance of OpenMP
%   testTracking                        - Test results of tracking against values previously computed
%   trackingdemo                        - Self-running tutorial
%   trackMultipleParticles              - Time the tracking of multiple particles
%   
%   AT/ATMAT/ATDEMOS/TRACKWITHIMPEDANCE
%   makeBBR                             - SI units;
%   testTrackBBR                        - Create a fast ring from ESRF lattice
%   
%   AT/ATMAT/ATGUI
%   intelem                             - Interactive element editor
%   intlat                              - Interactive AT lattice editor
%   
%   AT/ATMAT/ATMATCH
%   atApplyVariation                    - Applies a series of variations to variables as described in
%   atDisplayConstraintsChange          - This funciton evaluates the contraints defined in Constraints for lattice
%   atDisplayVariableChange             - This functions retrives variable Values for two rings to compare
%   atEvaluateConstraints               - This funciton evaluates the contraints defined in Constraints for lattice
%   atGetPenalty
%   atGetVariableValue
%   atlinconstraint                     - Generate a constraint on linear optics for atmatch
%   atmatch                             - Function [..
%   atVariableBuilder                   - AtVarBuilder   create a simple variable structure for use with atmatch
%   
%   AT/ATMAT/ATPHYSICS
%   atavedata                           - Average of optical functions on selected elements
%   findrespm                           - Computes the change in the closed orbit due to parameter perturbations
%   findspos                            - Returns longitudinal positions of accelerator lattice elements
%   PhysConstant                        - Physical constants
%   
%   AT/ATMAT/ATPHYSICS/COLLECTIVEEFFECTS
%   atbeam                              - Generates a particle distribution according to a sigma matrix
%   atsigma                             - Constructs a beam sigma matrix 2x2 4x4 or 6x6
%   
%   AT/ATMAT/ATPHYSICS/LINEAROPTICS
%   amat                                - Find A matrix from one turn map matrix T such that:
%   atdampingrates                      - Find tunes and damping rates from one map matrix with radiation
%   atlinopt                            - Performs 4D linear analysis of the COUPLED lattices
%   atlinopt2                           - Performs the linear analysis of UNCOUPLED lattices
%   atlinopt4                           - Performs the 4D linear analysis of COUPLED lattices
%   atlinopt6                           - Performs linear analysis of the lattice
%   beam22                              - Computes the beam matrix from the 1-turn transfer matrix
%   beam44                              - Computes the coupled beam matrices
%   find_betaoids                       - [H1 H2 H3]=find_betaoids(A)
%   find_etaoids                        - Given the normalizing matrix A, we compute the etaoids
%   find_inv_G                          - This function computes the invariants of a one turn map matrix
%   find_inv_G_fromA                    - This function computes the invariant matrices G1,G2,G3
%   findelemm44                         - Numerically finds the 4x4 transfer matrix of an element
%   findelemm66                         - Numerically finds the 6x6 transfer matrix of an element
%   findm44                             - Numerically finds the 4x4 transfer matrix of an accelerator lattice
%   findm66                             - Numerically finds the 6x6 transfer matrix of an accelerator lattice
%   get_dispersion_from_etaoids         - Computes dispersion functions (x,px,y,py) at refpts
%   jmat                                - Compute antisymmetric Matrix [O 1; -1 0]
%   linopt                              - Performs linear analysis of the COUPLED lattices
%   mkSRotationMatrix                   - (PSI) coordinate transformation matrix
%   plotbeta                            - Plots UNCOUPLED! beta-functions
%   r_analysis                          - Return the phase for A standardization
%   
%   AT/ATMAT/ATPHYSICS/LONGITUDINALDYNAMICS
%   atBunchLength                       - Bunch length due to the potential well effect
%   atRFacc                             - Computes RF acceptance of the ring
%   atsetcavity                         - ATSECAVITY Set the parameters of RF cavities
%   atSetCavityPhase                    - SETCAVITYPHASE     Set the TimeLag attribute of RF cavities
%   BunchLength                         - Bunch length due to the potential well effect
%   cavityoff                           - Turns cavities OFF
%   cavityon                            - Turns Cavities ON
%   mcf                                 - Momentum compaction factor
%   nus                                 - Computes synchrotron tune from RF parameters
%   phis                                - Phase = phis(U0MeV,VrfMV)
%   RFacc                               - Computes the RF acceptance with linear formula
%   
%   AT/ATMAT/ATPHYSICS/NAFFLIB
%   calcnaff                            - Computes NAFF decomposition for a phase space trajectory
%   naff_cc                             - Compile nafflibrary for Matlab
%   naff_example                        - Example to test naff within matlab
%   nafflib                             - MATLAB to NAFF library
%   
%   AT/ATMAT/ATPHYSICS/NONLINEARDYNAMICS
%   atnuampl                            - Computes tune shift with amplitude
%   computeRDT                          - Computes Hamiltonian resonance driving terms (RDTs)
%   tunespaceplot                       - Draws a tune diagram
%   
%   AT/ATMAT/ATPHYSICS/ORBIT
%   findorbit                           - Find the closed orbit
%   findorbit4                          - Finds closed orbit in the 4-d transverse phase
%   findorbit6                          - Finds closed orbit in the full 6-d phase space
%   findsyncorbit                       - Finds closed orbit, synchronous with the RF cavity
%   plotcod                             - Closed Orbit Distortion
%   xorbit_6                            - Private function used by findorbit6
%   xorbit_ct                           - Private function used by findsyncorbit
%   xorbit_dp                           - Private function used by findorbit4
%   
%   AT/ATMAT/ATPHYSICS/PARAMETERSUMMARYFUNCTIONS
%   atsummary                           - Print out various parameters of the current AT lattice
%   atx                                 - Computes and displays global information
%   RadIntegrals                        - Calcuate the contribution to the radiation integrals of a Wiggler
%   ringpara                            - Print out various parameters of the current AT lattice
%   twissline                           - Calculates linear optics functions for an UNCOUPLED transport line
%   twissring                           - Calculates linear optics functions for an UNCOUPLED ring
%   
%   AT/ATMAT/ATPHYSICS/RADIATION
%   atdisable_6d                        - Switches radiation and cavity off
%   atenable_6d                         - Switches RF and radiation on
%   atenergy                            - Gets the lattice energy
%   atgetU0                             - Computes Energy loss per turn in eV 
%   atradoff                            - Obsolete: switches RF and radiation off
%   atradon                             - Obsolete: switches RF and radiation on
%   atsetenergy                         - (ring,Energy) sets the Energy field in all
%   attapering                          - Scale magnet strengths
%   check_6d                            - Check the presence of longitudinal motion in a lattice
%   check_radiation                     - Obsolete: check the radiation state of a ring
%   DipoleRadiation                     - Compute the radiation integrals in dipoles
%   ElementRadiation                    - - Compute the radiation integrals in dipoles
%   ElossRadiation                      - Compute the radiation integrals in EnergyLoss elements
%   findelemraddiffm
%   findmpoleraddiffmatrix              - #ok<STOUT>
%   findthickmpoleraddiffm
%   findthinmpoleraddiffm
%   getclass_6d                         - Private. Guess class for 6d motion
%   ohmienvelope                        - Calculates equilibrium beam envelope in a
%   quantumDiff                         - Compute the radiation-diffusion matrix
%   radiationoff                        - Turns classical radiation  OFF
%   radiationon                         - Turns classical radiation  ON
%   thickmpoleraddiffm                  - FINDTHICKMPOLERADDIFFM
%   thinmpoleraddiffm                   - FINDTHINMPOLERADDIFFM
%   WigglerRadiation                    - Compute the radiation integrals in wigglers
%   
%   AT/ATMAT/ATPHYSICS/TOUSCHEKPIWINSKI
%   MomAperture_allRing                 - All Ring momentum aperture
%   momentum_aperture_at                - Function [deltamin, deltamax..
%   simpletestToucheckLT                - Simple example of use of toucheck lifetime formula:
%   TLT_IntPiw                          - Integral in Piwinski Formula for the Lifetime
%   TLT_IntPiw_k                        - Integral in Piwinski Formula for the Lifetime with u=tan^2(k)
%   TouschekPiwinskiLifeTime            - Function [Tl,contributionsTL]=TouschekPiwinskiLifeTime(ring,dpp,Ib,...)
%   
%   AT/ATMAT/ATPHYSICS/TUNEANDCHROMATICITY
%   findtune                            - Get the tune value from turn by turn positions
%   fitchrom2                           - Fits chromaticity  of THERING using 2 sextupole families
%   fittune2                            - Fits linear tunes of THERING using 2 quadrupole families
%   intfft                              - Calculates the tune from interpolated FFT of the trajectory
%   tunechrom                           - Computes linear tunes and chromaticities
%   
%   AT/ATMAT/ATPLOT
%   atbaseplot                          - Plots data generated by a user-supplied function
%   atplot                              - Plots optical functions
%   atplotsyn                           - Helper function for ATPLOT
%   atreforbit                          - Keep track of the nominal reference orbit through displaced elements
%   xplot                               - Private function used by atplot and atbaseplot
%   
%   AT/ATMAT/ATPLOT/PLOTFUNCTIONS
%   CurlyH                              - Function [H,Hv]=CurlyH(RING,dp,ind)
%   CurlyHlindata                       - Function [H,Hv]=CurlyHlindata(lindata)
%   plBeamSize                          - Plot H and V beam size
%   plClosedOrbit                       - Plots H and V 4x4 closed orbit
%   plCorrectorStrength                 - Plot PolynomB
%   plEmitContrib                       - Plot H/rho³ at every dipole
%   plenvelope                          - Plot beam envelope
%   plot_betabeat                       - Function plot_betabeat(THERING_ref,THERING_mod)
%   plot_trajectory                     - Plots particle trajectories
%   plotAperture                        - Plots x and y aperture
%   plotB0curlyh                        - Plot B and H
%   plotbetadisp                        - Function [s,plotdata]=plotbetadisp(ring,dpp,plotfun,varargin)
%   plotbetadispcurlyh                  - Plot beta, dispersion and H
%   plotERAperture                      - Plot RApertures EApertures
%   plotRDT
%   plotsqrtbetadispcurlyh              - Plot sqrt(beta), dispersion and H
%   plotWdispP                          - Plot W functions
%   plPolynomBComp                      - PlotBn coefficient with normalization
%   plPolynomBSxtOct                    - Plots Bn for sextupole and octupole magnets
%   plSigmaSigmap                       - Plots beam sizes and divergences
%   pltouschek                          - Plots Touschek lifetime contribution
%   plxi                                - #ok<INUSD>
%   
%   AT/ATMAT/ATTESTS
%   githubrun                           - Private. This script runs the AT test suite in a GitHub action. There is
%   githubsetup                         - Private. Setup Matlab for AT tests in GitHib Actions
%   pytests                             - Shared setup for the entire test class
%   
%   AT/ATMAT/ATTRACK
%   atpass                              - #ok<STOUT>
%   linepass                            - Tracks particles through each element of the cell array LINE
%   ringpass                            - Tracks particles through each element of the cell array RING
%   
%   AT/ATMAT/ATUTILS
%   atoptions                           - Definition of default parameters
%   frequency_control                   - Private. Handle off-momentum for 6D lattice
%   getargs                             - Process positional arguments from the input arguments
%   getdparg                            - Handle positional dp arguments
%   getenvopt                           - (NAME, DEFAULTVALUE)
%   getflag                             - Check the presence of a flag in an argument list
%   getoption                           - Extract a keyword argument from an argument list
%   opticsoptions                       - (private) extract arguments for atlinopt
%   parseargs                           - Check and expands optional argument lists
%   setoption                           - Set AT preference values
%   wrapper6d                           - Private. Handle off-momentum for 6D lattice
%   
%   AT/ATMAT/LATTICE
%   at2py                               - ELSTR=AT2PY(ELEM) convert AT element tp pyat
%   at2str                              - Makes the string representation of an AT element
%   ataddmpolecomppoly                  - Adds a multipole component to an existing polynomial,
%   ataddmpoleerrors                    - Ataddrandmpole adds a random multipole component to all elements of type
%   atCheckRingProperties               - Get the ring properties if existing
%   atdivelem                           - Divide an element into pieces
%   atelem                              - Makes a new AT element structure from another element,
%   atfastring                          - Generate simplified AT structures
%   atfitchrom                          - Fit chromaticites by scaling 2 sextupole families
%   atfittune                           - Fit linear tunes by scaling 2 quadrupole families
%   atgetcells                          - Performs a search on MATLAB cell arrays of structures
%   atgetfieldvalues                    - Retrieves the field values AT cell array of elements
%   atGetRingProperties                 - Get the ring properties
%   atguessclass                        - Tries to determine the class of an element
%   atindex                             - Extracts the information about element families and
%   atinsertelems                       - Insert elements at given locations in a line
%   atloadfielderrs                     - Will load a field error structure into a ring
%   atloadlattice                       - Load a lattice from a file
%   atlocateparam                       - Private function. Locate the RingParam element
%   atmaincavities                      - Get the fundamental mode cavities
%   atmakefielderrstruct                - MAKERNDFIELDERRS will create a field error data structure
%   atparamscan                         - Private function. Updates the RingParam element
%   atparticle                          - Particle definition for AT
%   atreduce                            - Remove useless elements from an AT structure
%   atrotatelattice                     - Circularly shift the lattice elements
%   atsbreak                            - Insert markers at given s positions in a lattice
%   atsetfieldvalues                    - Sets the field values of MATLAB cell array of structures
%   atsetRFCavity                       - Set the RF Cavity with the passmethod RFCavityPass
%   atSetRingProperties                 - Add or modify properties of the lattice
%   atsetshift                          - Sets the misalignment vectors
%   atsettilt                           - Sets the entrance and exit rotation matrices
%   atshiftelem                         - Set new displacement parameters
%   atsimplering                        - Creates a "simple ring"
%   atsplitelem                         - Creates a line by inserting one or more elements into a base element
%   attiltelem                          - Sets new rotation parameters
%   atwritem                            - Creates a .m file to store an AT structure
%   atwritepy                           - Creates pyAT lattice from a Matlab lattice
%   buildlat                            - Places elements from FAMLIST into cell array THERING
%   combinebypassmethod                 - Combines adjacent elements that have the same specified pass method
%   combinelinear45                     - Combines adjacent  elements that use 4-by-5 PassMethods
%   entrancefields                      - () Return the list of field names affecting the element entrance
%   exitfields                          - () Return the list of field names affecting the element exit
%   findcells                           - Performs a search on MATLAB cell arrays of structures
%   findtags                            - Looks for string matches in 'Tag' field of AT lattice elements
%   getcellstruct                       - Retrieves the field values MATLAB cell array of structures
%   insertelem0                         - - quick and dirty:
%   insertindrift                       - Inserts one or more elements into a drift element
%   isatelem                            - Tests if an input argument is a valid AT element
%   isatlattice                         - Tests if an input argument is a valid AT lattice
%   mergedrift                          - Removes a lattice element and merges the two adjacent drift spaces
%   mvelem                              - Move an element
%   mvfield                             - Move fields from one structure to another
%   rmelem0                             - Removes elements of length 0 from the accelerator lattice
%   setcellstruct                       - Sets the field values of MATLAB cell array of structures
%   setshift                            - Sets the misalignment vectors T1, T2 for elements
%   settags                             - Sets the 'Tag' field in AT lattice elements
%   settilt                             - Sets the entrance and exit misalignment matrixes
%   splitdrift                          - Inserts an element into a drift space
%   symplectify                         - Makes a matrix more symplectic
%   
%   AT/ATMAT/LATTICE/CONVERTERS
%   readmad                             - Reads the file output of MAD commands
%   
%   AT/ATMAT/LATTICE/CONVERTERS/AT2ELEGANT
%   AT_2_Elegant                        - This functions converts the AT lattice AT_ring in elegant form
%   
%   AT/ATMAT/LATTICE/CONVERTERS/AT2G4BL
%   ATtoG4BL                            - Function [outtext]=ATtoG4BL(P_0,particle,folder)
%   
%   AT/ATMAT/LATTICE/CONVERTERS/AT2MAD8
%   AT_2_mad8                           - Function [elelat,def,lines]=AT_2_mad8(AT_ring,linename)
%   
%   AT/ATMAT/LATTICE/CONVERTERS/AT2MADX
%   AT_2_madX                           - Function [elelat,defs,lines]=AT_2_madX(AT_ring,linename)
%   
%   AT/ATMAT/LATTICE/CONVERTERS/AT2OPA
%   AT_2_OPA                            - Function AT_2_OPA(AT_ring,linename)
%   
%   AT/ATMAT/LATTICE/CONVERTERS/ELEGANT2AT
%   ele2at_run_me                       - Test_elegant_converter
%   elegant2at                          - Function elegant2at(elegantlattice,E0,outfilename)
%   ParseAtributesELEGANT_2_AT          - Determines atribute and sets field in sxs{i} structure AT
%   
%   AT/ATMAT/LATTICE/CONVERTERS/MAD82MADX
%   mad8TOmadx                          - Converts mad8 sequence files to madX
%   
%   AT/ATMAT/LATTICE/CONVERTERS/MADX2AT
%   atfrommadx                          - Function atfrommadx(seqfilemadX,E0,outfilename)
%   buildATLattice                      - Given a list (cell array) of elements with specified field Spos (center of element (madx default)) in a
%   ParseAtributesMADX_2_AT             - Determines atribute and sets field in sxs{i} structure AT
%   reshapeToCellArray                  - If CEL_CEL is a cell array of structures and cell arrays it converts it a
%   
%   AT/ATMAT/LATTICE/CONVERTERS/MADX2AT/EXAMPLES
%   convertMADXtoATExample              - Simple lattice test (uncomment to run)
%   
%   AT/ATMAT/LATTICE/CONVERTERS/MADX2G4BL
%   madx2g4bl                           - Function [outtext]=madx2g4bl(P_0,particle,folder)
%   
%   AT/ATMAT/LATTICE/ELEMENT_CREATION
%   ataperture                          - Creates a aperture element
%   atbaselem                           - Create an AT element structure + various checks
%   atcorrector                         - Creates a drift space element with class 'Corrector'
%   atdampMatElem                       - Creates an element that applies the global damping matrix
%   atdrift                             - Creates a drift space element with Class 'Drift'
%   atenergyloss                        - Creates an energy loss element
%   atidtable                           - Creates an ID element
%   atinsertiondevicekickmap            - Creates an insertion device kick-map element
%   atM66                               - Create an element applying an arbitrary 6x6 transfer matrix
%   atM66Tijk                           - ATM66(FAMNAME,M66,Tijk,PASSMETHOD)
%   atmarker                            - Creates a marker space element
%   atmonitor                           - Creates a Beam Position Monitor element with Class 'Monitor'
%   atmultipole                         - Creates a multipole element
%   atquadrupole                        - Creates a quadrupole element with Class 'Quadrupole'
%   atQuantDiff                         - Creates a quantum diffusion element
%   atrbend                             - Creates a rectangular bending magnet element with class 'Bend'
%   atrbendtune                         - Set X0ref and RefDZ for rectangular bending magnets
%   atrfcavity                          - Creates an rfcavity element with Class 'RFCavity'
%   atringparam                         - Creates a RingParameter Element which should go at the beginning of the ring
%   atsbend                             - Creates a sector bending magnet element with class 'Bend'
%   atsextupole                         - Creates a sextupole element with class 'Sextupole'
%   atSimpleQuantDiff                   - SimpleQuantDiff creates a simple quantum difusion element
%   atskewquad                          - Creates a skew quadrupole element with Class 'Multipole'
%   atsolenoid                          - Creates a new solenoid element with Class 'Solenoid'
%   atthinmultipole                     - Creates a thin multipole element
%   atvariablemultipole                 - Creates a variable thin multipole element
%   atwiggler                           - Creates a wiggler
%   
%   AT/ATMAT/LATTICE/ELEMENT_CREATION/PRIVATE
%   decodeatargs                        - Separates arguments and resources
%   
%   AT/ATMAT/LATTICE/PARAMGROUP
%   atparamgroup                        - PARAMETER GROUP in AT is a general way
%   mkparamgroup                        - Simplifies creation of AT parameter groups
%   restoreparamgroup                   - Restores the values of multiple physical
%   saveparamgroup                      - Saves the values of multiple physical
%   setparamgroup                       - Modifies a group of parameters
%   
%   AT/ATMAT/LATTICE/SURVEY
%   atgeometry                          - Computes the 2-D position of all elements (no vertical bend)
%   atgeometry3                         - Computes the 3-D position of all elements
%   
%   AT/ATMAT/PUBTOOLS
%   atdynap                             - Compute the dynamic aperture
%   atmomap                             - Find momentum aperture at start of ring
%   atsurvey2spos                       - Returns closest lattics s coordinates to xycoord points
%   atundulator                         - Define undulator model
%   atvalue                             - Extract array from lindata structure
%   calc_dppAperture                    - Calculate the momentum aperture at each location of the ring due to
%   calc_Touschek                       - TauT = calc_Touschek(THERING, Ib)
%   calc_TouschekPM                     - TauT = calc_TouschekPM(TD,dppPM,Trf,Ib,U0,coupling, sigE, emit_x)
%   freqsearch                          - =========================================================================
%   nlchromplot                         - Example:   nlchromplot(esrf,-.04,.04,30,16,1)
%   
%   AT/ATMAT/PUBTOOLS/APERTURE
%   SetPhysicalAperture                 - Ringapert=SetPhysicalAperture(ring,apertureX,apertureY)
%   
%   AT/ATMAT/PUBTOOLS/CREATE_ELEMS
%   atidtable_dat                       - Atidtable(FamName, Nslice, filename, Energy, method)
%   
%   AT/ATMAT/PUBTOOLS/DISTANCE2CURVE
%   distance2curve                      - Gets the minimum distance from a point to a general curvilinear n-dimensional arc
%   
%   AT/ATMAT/PUBTOOLS/HAISSINSKI
%   blength
%   fitgaussian                         - GAUSSIAN_PARAM FITERR GAUSSFIT SIGERROR]= FITGAUSSIAN(DATA,[property_value_pair]);
%   plothaissinski                      - [z lambda sigma mu] = HASSINSKYFIT(SIGMA0, R, L)
%   Qval                                - Gives the unitless Q parameter needed to compute the
%   
%   AT/ATMAT/PUBTOOLS/LATTICE_TOOLS
%   atreadbeta                          - Reads a BETA file
%   atsetglobval                        - Creates the global variable GLOBVAL and adds Energy
%   scalesext                           - Newring=scalesext(ring,sextfam,scale)
%   setsext                             - Newring=setsext(ring,fam,val);
%   setsextall                          - Newring=setsextall(ring,fam,val);
%   sext_sens_scan                      - Esrf(:)
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION
%   qemsvd_mod                          - Function dq=qemsvd_mod(a,b,neig,plot)
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/CHROMATICITY
%   atmatchchromdelta                   - Function arcchrom0=atmatchchromdelta(arc,c,sxtfams)
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/CORRECTION_CHAIN
%   CorrectionChain                     - Corrected lattice
%   DisplayCorrectionEffect             - [d0,de,dc]=DisplayCorrectionEffect(..
%   testcorrectionchain                 - Test errors and correction functions
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/DISPERSION
%   atcorrectdispersion                 - Function [..
%   testdispersioncorrection            - Test errors and correction functions
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/DISPERSIONFREESTEERING
%   atdispersionfreesteering            - Function [..
%   testdispersionfreesteering          - Test errors and correction functions
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/ORBIT
%   atcorrectorbit
%   testorbitbump                       - Test errors and correction functions
%   testorbitcorrection                 - Test errors and correction functions
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/ORBITBUMPS/MATCHING
%   BumpAtBPM                           - Function roff=BumpAtBPM(..
%   BumpAtBPM4D                         - Function roff=BumpAtBPM(..
%   testorbitbump                       - Test matching orbit bump
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/RDT
%   atavedata_mod                       - ATAVEDATA       Average of optical functions on selected elements
%   atQuadRDTdispersioncorrection       - - Make dispersion correction based on RDTs
%   atRDTdispersioncorrection           - Makes dispersion correction based on RDTs
%   atRDTdispersionmeasuredcorrection   - Makes correction of dispersion based on
%   atSkewRDTdispersioncorrection       - Function [..
%   EquivalentGradientsFromAlignments6D - Estimated normal quad gradients from sext offsets
%   qemrdtresp_mod                      - QEMRDTRESP  compute resonance driving terms at BPM locations
%   semrdtresp_mod                      - SEMRDT compute resonance driving terms at BPM locations
%   testRDTdispersionfreesteering       - Test errors and correction functions
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/RESPONSE_MATRIX
%   findrespmat                         - FINDRESPM_mod computes the change in the closed orbit due to parameter perturbations
%   getresponsematrices                 - 1 AT lattice
%   gettunechromatlinopt                - Gets tunes and chromaticities from atlinopt
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/RFCAVITY
%   atRFcorrection                      - Function [..
%   atsetRFCavityErr                    - ATSETRFCAVITY sets the RF Cavity with the passmethod RFCavityPass
%   testsetRFCavityErr                  - Test errors and correction functions
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/TRAJECTORY
%   atfirstturntrajectory               - Makes first turn correction
%   MatchLast2CorForFirstBPM            - Takes the last two correctors to match the orbit and angle trajectory at
%   Scan2x2DinCOD                       - [bestinputcoord]=ScanPosAngle(..
%   testorbitaftertrajectory            - Test errors and correction functions
%   testtrajectorycorrection            - Test errors and correction functions
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/TUNE
%   atmatchtunedelta                    - Function arcchrom0=atmatchtunedelta(arc,c,quadfams)
%   fittunedelta2fam                    - Rerr=fittunedelta2fam(rerr,r0)
%   testfittunedelta2fam                - Test errors and correction functions
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS
%   AssignFieldErr                      - Function r=AssignFieldErr(r,refpos,N,rho,BNn,ANn)
%   atset_s_shift                       - Implements DS longitudinal position drift
%   atsetbpmerr                         - Sets the misalignment vectors
%   atsettiltdipole                     - Sets the entrance and exit rotation matrices
%   bpm_matrices                        - Generate transformation matrices for BPM readings
%   bpm_process                         - Compute BPM readings from the closed orbit
%   finddispersion6Err                  - Gets 6D dispersion with bpm reading errors
%   findorbit4Err                       - Gets 4x4 closed orbit with BPM errors
%   findorbit6Err                       - Findorbit6 with bpm reading errors
%   findtrajectory6Err                  - [t    6xNbpm array of  trajectory
%   setANYshift                         - Adds to the existing shift errors additional D
%   setFieldIntegralError               - Function rerr=setFieldIntegralError(r0,rerr,indx,order,Nsigma,sigmaperc)
%   setGirderError                      - Rerr=setGirderError(r,pert,errval,mag_group)
%   SetLargeErrorList                   - Sets given error list
%   setTiltAbout                        - Sets tilt errors
%   setTiltGirderAbout                  - Set Tilt error on a magnet
%   setXshift                           - Set horizontal shifts for a element list
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/BPMERRORS
%   testerr                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/DELTAS
%   testerr                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/DXDY
%   testerr                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/ERRORDISPLAYFUNCTIONS
%   GetMisalignments                    - This function retrives 3 vectors, for x and y misalignments and tilts
%   pltmisalignments                    - #ok<INUSD>
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/ERRORSMANIPULATION
%   GetExistingErrors                   - This function retrives 6 vectors, for x, y,s misalignments,
%   getMagGroupsFromGirderIndex         - Gets magnets on a girder
%   getMagGroupsFromMagNum              - - Gets magnet from a Magnet group
%   setBpmOffsetOnDipoleRef             - Set bpm on curve defined by dipole misalignments
%   SetExistingError                    - Function SetExistingError(rerr,magindex,X0,Y0,S0,T0,R0,P0,bpm0)
%   SumErrors                           - Rsum=SumErrors(r1,r2,magindex)
%   ThetaPhiGirder                      - Rtp=ThetaPhiGirder(r,mag_gr)
%   UniformGirderErrors                 - Function ring=UniformGirderErrors(ring)
%   UniformMagGroupsErrors              - Function ring=UniformMagGroupsErrors(ring)
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/FIELDINTEGRAL
%   testerr                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/GIRDERS
%   testerr                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/LARGEERRLIST
%   testerr                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/LARGEWAVELIST
%   testcor                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/MULTIPOLES
%   testerr                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/RANDOM
%   atsetrandomerrors                   - Function rerr=atsetrandomerrors(..
%   seterrorrand                        - Nominal lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/SURVEY
%   SetESRFAlgeAlignmentError           - Function SetESRFAlgeAlignmentError(..
%   testerr                             - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/TILT
%   compRotCorVsDip                     - Load lattice
%   compRotCorVsDipQuad                 - Load lattice
%   testerr                             - Load lattice
%   testerrRotCorrector                 - Load lattice
%   
%   AT/ATMAT/PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/WAVE
%   atsetwaveerrors                     - Function rerr=atsetwaveerrors(..
%   seterrorwave                        - Nominal lattice
%   
%   AT/ATMAT/PUBTOOLS/LOCAL_LATTICE_PARAMS
%   atmakeXYProjectionEllipse           - Gives points to plot the contour ellipses
%   machine_at                          - Machine AT will return the optics of the lattice. Essentially takes what
%   plotContours                        - Plots contours
%   
%   AT/ATMAT/PUBTOOLS/VACUUMLIFETIME
%   VacLifetimeResidualGass             - Coloumb scattering and residual gas Bremsstrahlung Cross sections are
%   
%   AT/ATOCTAVE
%   atmexall                            - Build all AT platform dependent mex-files from C-sources
%   atparticle                          - Particle definition for AT
%   bootstrap                           - Prepare Octave for AT
%   endsWith                            - True if text ends with pattern
%   isstring                            - Determine whether input is string array
%   octaveVersion                       - Get major and minor version numbers
%   startsWith                          - True if text starts with pattern
%   
%   AT/BUILD/LIB.MACOSX-13-X86_64-CPYTHON-39/MACHINE_DATA
%   macosx-13-x86_64-cpython-39/machine_data.australian_synchrotron
%   macosx-13-x86_64-cpython-39/machine_data.esrf
%   macosx-13-x86_64-cpython-39/machine_data.soleil
%   macosx-13-x86_64-cpython-39/machine_data.thomx
%   
%   AT/DEVELOPER/MATLAB
%   atchapters                          - Describe the User Guide chapters
%   atclearmex                          - Remove all AT mex-files
%   atrelease                           - RELEASE    build, test and package a the AT package
%   gen_help                            - Build the "help" infrastructure
%   gen_list                            - Display a list of still non-documented AT functions
%   gen_toc                             - Build the HTML files used by the Matlab help browser
%   h1_line                             - Get the H1 line for a file
%   howtochapters                       - Describe the How to chapters
%   setversion                          - Set the version of AT
%   
%   AT/DEVELOPER/MATLAB/M
%   atelemcreate                        - Element creation
%   atelemfuncs                         - Element manipulation
%   atlatticefuncs                      - Lattice manipulation
%   atlinearoptics                      - Linear optics
%   atloadsave                          - Loading and Saving lattices
%   atphysics                           - Physics
%   atradiation                         - Radiation
%   atsummary                           - Parameter summary
%   howtosummary                        - How to…
%   ugsummary                           - AT User Guide
%   
%   AT/DOC
%   talk
%   
%   AT/MACHINE_DATA
%   australian_synchrotron
%   dba                                 - Create dba lattice
%   esrf
%   FODO                                - P1Dr=atdrift('Dr',0.1);
%   soleil                              - Loads SOLEIL lattice
%   sp3v81f                             - All the dipole and quadrupole lengths are effective lengths
%   spear2                              - Example lattice definition file
%   spear2rad                           - Example lattice definition file with CAVITY and CLASSICAL radiation
%   spear2resp                          - Example SPEAR2 lattice with orbit correctors and BPMS
%   spear3                              - Load the SPEAR3 lattice structure
%   thomx
%   
%   AT/PYAT/MACHINE_DATA
%   australian_synchrotron
%   esrf
%   soleil
%   thomx
%   
%   AT/PYAT/TEST_MATLAB
%   pyproxy                             - Convert structures to a form accessible in python
%   
%   AT/UTILS/MPI_SWEEP
%   mpi_sweep_octave
%   mpi_sweep_octave_example            - D%d", i, j);
%    
%   This file was generated by updateContents.m on 06 Nov 2023 at 14:58:15.



