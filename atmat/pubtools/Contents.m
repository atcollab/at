% PUBTOOLS
% See also 
% 
%   Contents file for PUBTOOLS and its subfolders.
%   
%   PUBTOOLS
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
%   PUBTOOLS/APERTURE
%   SetPhysicalAperture                 - Ringapert=SetPhysicalAperture(ring,apertureX,apertureY)
%   
%   PUBTOOLS/CREATE_ELEMS
%   atidtable_dat                       - Atidtable(FamName, Nslice, filename, Energy, method)
%   
%   PUBTOOLS/DISTANCE2CURVE
%   distance2curve                      - Gets the minimum distance from a point to a general curvilinear n-dimensional arc
%   
%   PUBTOOLS/HAISSINSKI
%   blength
%   fitgaussian                         - GAUSSIAN_PARAM FITERR GAUSSFIT SIGERROR]= FITGAUSSIAN(DATA,[property_value_pair]);
%   plothaissinski                      - [z lambda sigma mu] = HASSINSKYFIT(SIGMA0, R, L)
%   Qval                                - Gives the unitless Q parameter needed to compute the
%   
%   PUBTOOLS/LATTICE_TOOLS
%   atreadbeta                          - Reads a BETA file
%   atsetglobval                        - Creates the global variable GLOBVAL and adds Energy
%   scalesext                           - Newring=scalesext(ring,sextfam,scale)
%   setsext                             - Newring=setsext(ring,fam,val);
%   setsextall                          - Newring=setsextall(ring,fam,val);
%   sext_sens_scan                      - Esrf(:)
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION
%   qemsvd_mod                          - Function dq=qemsvd_mod(a,b,neig,plot)
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/CHROMATICITY
%   atmatchchromdelta                   - Function arcchrom0=atmatchchromdelta(arc,c,sxtfams)
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/CORRECTION_CHAIN
%   CorrectionChain                     - Corrected lattice
%   DisplayCorrectionEffect             - [d0,de,dc]=DisplayCorrectionEffect(..
%   testcorrectionchain                 - Test errors and correction functions
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/DISPERSION
%   atcorrectdispersion                 - Function [..
%   testdispersioncorrection            - Test errors and correction functions
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/DISPERSIONFREESTEERING
%   atdispersionfreesteering            - Function [..
%   testdispersionfreesteering          - Test errors and correction functions
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/ORBIT
%   atcorrectorbit
%   testorbitbump                       - Test errors and correction functions
%   testorbitcorrection                 - Test errors and correction functions
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/ORBITBUMPS/MATCHING
%   BumpAtBPM                           - Function roff=BumpAtBPM(..
%   BumpAtBPM4D                         - Function roff=BumpAtBPM(..
%   testorbitbump                       - Test matching orbit bump
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/RDT
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
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/RESPONSE_MATRIX
%   findrespmat                         - FINDRESPM_mod computes the change in the closed orbit due to parameter perturbations
%   getresponsematrices                 - 1 AT lattice
%   gettunechromatlinopt                - Gets tunes and chromaticities from atlinopt
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/RFCAVITY
%   atRFcorrection                      - Function [..
%   atsetRFCavityErr                    - ATSETRFCAVITY sets the RF Cavity with the passmethod RFCavityPass
%   testsetRFCavityErr                  - Test errors and correction functions
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/TRAJECTORY
%   atfirstturntrajectory               - Makes first turn correction
%   MatchLast2CorForFirstBPM            - Takes the last two correctors to match the orbit and angle trajectory at
%   Scan2x2DinCOD                       - [bestinputcoord]=ScanPosAngle(..
%   testorbitaftertrajectory            - Test errors and correction functions
%   testtrajectorycorrection            - Test errors and correction functions
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/CORRECTION/TUNE
%   atmatchtunedelta                    - Function arcchrom0=atmatchtunedelta(arc,c,quadfams)
%   fittunedelta2fam                    - Rerr=fittunedelta2fam(rerr,r0)
%   testfittunedelta2fam                - Test errors and correction functions
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS
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
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/BPMERRORS
%   testerr                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/DELTAS
%   testerr                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/DXDY
%   testerr                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/ERRORDISPLAYFUNCTIONS
%   GetMisalignments                    - This function retrives 3 vectors, for x and y misalignments and tilts
%   pltmisalignments                    - #ok<INUSD>
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/ERRORSMANIPULATION
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
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/FIELDINTEGRAL
%   testerr                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/GIRDERS
%   testerr                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/LARGEERRLIST
%   testerr                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/LARGEWAVELIST
%   testcor                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/MULTIPOLES
%   testerr                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/RANDOM
%   atsetrandomerrors                   - Function rerr=atsetrandomerrors(..
%   seterrorrand                        - Nominal lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/SURVEY
%   SetESRFAlgeAlignmentError           - Function SetESRFAlgeAlignmentError(..
%   testerr                             - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/TILT
%   compRotCorVsDip                     - Load lattice
%   compRotCorVsDipQuad                 - Load lattice
%   testerr                             - Load lattice
%   testerrRotCorrector                 - Load lattice
%   
%   PUBTOOLS/LATTICETUNINGFUNCTIONS/ERRORS/WAVE
%   atsetwaveerrors                     - Function rerr=atsetwaveerrors(..
%   seterrorwave                        - Nominal lattice
%   
%   PUBTOOLS/LOCAL_LATTICE_PARAMS
%   atmakeXYProjectionEllipse           - Gives points to plot the contour ellipses
%   machine_at                          - Machine AT will return the optics of the lattice. Essentially takes what
%   plotContours                        - Plots contours
%   
%   PUBTOOLS/VACUUMLIFETIME
%   VacLifetimeResidualGass             - Coloumb scattering and residual gas Bremsstrahlung Cross sections are
%    
%   This file was generated by updateContents.m on 06 Nov 2023 at 14:58:13.
