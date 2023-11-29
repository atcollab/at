% ATPHYSICS
% See also 
% 
%   Contents file for ATPHYSICS and its subfolders.
%   
%   ATPHYSICS
%   atavedata                   - Average of optical functions on selected elements
%   findrespm                   - Computes the change in the closed orbit due to parameter perturbations
%   findspos                    - Returns longitudinal positions of accelerator lattice elements
%   PhysConstant                - Physical constants
%   
%   ATPHYSICS/COLLECTIVEEFFECTS
%   atbeam                      - Generates a particle distribution according to a sigma matrix
%   atsigma                     - Constructs a beam sigma matrix 2x2 4x4 or 6x6
%   
%   ATPHYSICS/LINEAROPTICS
%   amat                        - Find A matrix from one turn map matrix T such that:
%   atdampingrates              - Find tunes and damping rates from one map matrix with radiation
%   atlinopt                    - Performs 4D linear analysis of the COUPLED lattices
%   atlinopt2                   - Performs the linear analysis of UNCOUPLED lattices
%   atlinopt4                   - Performs the 4D linear analysis of COUPLED lattices
%   atlinopt6                   - Performs linear analysis of the lattice
%   beam22                      - Computes the beam matrix from the 1-turn transfer matrix
%   beam44                      - Computes the coupled beam matrices
%   find_betaoids               - [H1 H2 H3]=find_betaoids(A)
%   find_etaoids                - Given the normalizing matrix A, we compute the etaoids
%   find_inv_G                  - This function computes the invariants of a one turn map matrix
%   find_inv_G_fromA            - This function computes the invariant matrices G1,G2,G3
%   findelemm44                 - Numerically finds the 4x4 transfer matrix of an element
%   findelemm66                 - Numerically finds the 6x6 transfer matrix of an element
%   findm44                     - Numerically finds the 4x4 transfer matrix of an accelerator lattice
%   findm66                     - Numerically finds the 6x6 transfer matrix of an accelerator lattice
%   get_dispersion_from_etaoids - Computes dispersion functions (x,px,y,py) at refpts
%   jmat                        - Compute antisymmetric Matrix [O 1; -1 0]
%   linopt                      - Performs linear analysis of the COUPLED lattices
%   mkSRotationMatrix           - (PSI) coordinate transformation matrix
%   plotbeta                    - Plots UNCOUPLED! beta-functions
%   r_analysis                  - Return the phase for A standardization
%   
%   ATPHYSICS/LONGITUDINALDYNAMICS
%   atBunchLength               - Bunch length due to the potential well effect
%   atRFacc                     - Computes RF acceptance of the ring
%   atsetcavity                 - ATSECAVITY Set the parameters of RF cavities
%   atSetCavityPhase            - SETCAVITYPHASE     Set the TimeLag attribute of RF cavities
%   BunchLength                 - Bunch length due to the potential well effect
%   cavityoff                   - Turns cavities OFF
%   cavityon                    - Turns Cavities ON
%   mcf                         - Momentum compaction factor
%   nus                         - Computes synchrotron tune from RF parameters
%   phis                        - Phase = phis(U0MeV,VrfMV)
%   RFacc                       - Computes the RF acceptance with linear formula
%   
%   ATPHYSICS/NAFFLIB
%   calcnaff                    - Computes NAFF decomposition for a phase space trajectory
%   naff_cc                     - Compile nafflibrary for Matlab
%   naff_example                - Example to test naff within matlab
%   nafflib                     - MATLAB to NAFF library
%   
%   ATPHYSICS/NONLINEARDYNAMICS
%   atnuampl                    - Computes tune shift with amplitude
%   computeRDT                  - Computes Hamiltonian resonance driving terms (RDTs)
%   tunespaceplot               - Draws a tune diagram
%   
%   ATPHYSICS/ORBIT
%   findorbit                   - Find the closed orbit
%   findorbit4                  - Finds closed orbit in the 4-d transverse phase
%   findorbit6                  - Finds closed orbit in the full 6-d phase space
%   findsyncorbit               - Finds closed orbit, synchronous with the RF cavity
%   plotcod                     - Closed Orbit Distortion
%   xorbit_6                    - Private function used by findorbit6
%   xorbit_ct                   - Private function used by findsyncorbit
%   xorbit_dp                   - Private function used by findorbit4
%   
%   ATPHYSICS/PARAMETERSUMMARYFUNCTIONS
%   atsummary                   - Print out various parameters of the current AT lattice
%   atx                         - Computes and displays global information
%   RadIntegrals                - Calcuate the contribution to the radiation integrals of a Wiggler
%   ringpara                    - Print out various parameters of the current AT lattice
%   twissline                   - Calculates linear optics functions for an UNCOUPLED transport line
%   twissring                   - Calculates linear optics functions for an UNCOUPLED ring
%   
%   ATPHYSICS/RADIATION
%   atdisable_6d                - Switches radiation and cavity off
%   atenable_6d                 - Switches RF and radiation on
%   atenergy                    - Gets the lattice energy
%   atgetU0                     - Computes Energy loss per turn in eV 
%   atradoff                    - Obsolete: switches RF and radiation off
%   atradon                     - Obsolete: switches RF and radiation on
%   atsetenergy                 - (ring,Energy) sets the Energy field in all
%   attapering                  - Scale magnet strengths
%   check_6d                    - Check the presence of longitudinal motion in a lattice
%   check_radiation             - Obsolete: check the radiation state of a ring
%   DipoleRadiation             - Compute the radiation integrals in dipoles
%   ElementRadiation            - - Compute the radiation integrals in dipoles
%   ElossRadiation              - Compute the radiation integrals in EnergyLoss elements
%   findelemraddiffm
%   findmpoleraddiffmatrix      - #ok<STOUT>
%   findthickmpoleraddiffm
%   findthinmpoleraddiffm
%   getclass_6d                 - Private. Guess class for 6d motion
%   ohmienvelope                - Calculates equilibrium beam envelope in a
%   quantumDiff                 - Compute the radiation-diffusion matrix
%   radiationoff                - Turns classical radiation  OFF
%   radiationon                 - Turns classical radiation  ON
%   thickmpoleraddiffm          - FINDTHICKMPOLERADDIFFM
%   thinmpoleraddiffm           - FINDTHINMPOLERADDIFFM
%   WigglerRadiation            - Compute the radiation integrals in wigglers
%   
%   ATPHYSICS/TOUSCHEKPIWINSKI
%   MomAperture_allRing         - All Ring momentum aperture
%   momentum_aperture_at        - Function [deltamin, deltamax..
%   simpletestToucheckLT        - Simple example of use of toucheck lifetime formula:
%   TLT_IntPiw                  - Integral in Piwinski Formula for the Lifetime
%   TLT_IntPiw_k                - Integral in Piwinski Formula for the Lifetime with u=tan^2(k)
%   TouschekPiwinskiLifeTime    - Function [Tl,contributionsTL]=TouschekPiwinskiLifeTime(ring,dpp,Ib,...)
%   
%   ATPHYSICS/TUNEANDCHROMATICITY
%   findtune                    - Get the tune value from turn by turn positions
%   fitchrom2                   - Fits chromaticity  of THERING using 2 sextupole families
%   fittune2                    - Fits linear tunes of THERING using 2 quadrupole families
%   intfft                      - Calculates the tune from interpolated FFT of the trajectory
%   tunechrom                   - Computes linear tunes and chromaticities
%    
%   This file was generated by updateContents.m on 06 Nov 2023 at 14:58:12.
