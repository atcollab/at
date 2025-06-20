.. _parametersummaryfunctions_module:

ParameterSummaryFunctions
=========================

.. rubric:: Functions


.. list-table::

   * - :func:`atsummary`
     - ATSUMMARY  Print out various parameters of the current AT lattice
   * - :func:`twissring`
     - TWISSRING calculates linear optics functions for an UNCOUPLED ring
   * - :func:`RadIntegrals`
     - calcuate the contribution to the radiation integrals of a Wiggler.
   * - :func:`atx`
     - ATX Computes and displays global information
   * - :func:`twissline`
     - TWISSLINE calculates linear optics functions for an UNCOUPLED transport line
   * - :func:`ringpara`
     - RINGPARA  Print out various parameters of the current AT lattice

.. py:function:: atsummary

   | 
   | ATSUMMARY()                Print parameters of the global variable THERING
   | ATSUMMARY(RING)            Print parameters of RING. RING may be 4D or 6D
   | 
   | ATSUMMARY(...,'dp',DP)     Specify off-momentum. For off-momentum lattices,
   |    the equilibrium emittance and energy spread will not be computed
   |    because the quadrupole contribution is missing.
   | 
   | ATSUMMARY(...,'dct',DCT)   Specify the path lengthening
   | 
   | ATSUMMARY(...,'df',DF)     Specify the RF frequency deviation
   | 
   | RINGPARAMS=ATSUMMARY(...)  Return the results in a structure instead of
   |    printing them
   | 
   | RINGPARAMS=ATSUMMARY(...,'Display')  Restore the printout of parameters
   | 
   |   The parameters that come after the Synchrotron Integrals are
   |   parameters that depend on the Integrals themselves. The equations to
   |   calculate them were taken from [1].
   | 
   |   [1] Alexander Wu Chao and Maury Tigner, Handbook of Accelerator Physics
   |   and Engineering (World Scientific, Singapore, 1998), pp. 183-187. (or
   |   187-190 in ed. 2)
   | 
   | See also ATX RINGPARA

.. py:function:: twissring

   |  
   |  [TwissData, tune]  = TWISSRING(LATTICE,DP) calculates twiss parameters
   |     and closed orbit coordinates at the RING entrance assuming
   |     constant energy deviation DP.
   | 
   |  [TwissData, tune]  = TWISSRING(LATTICE,DP,REFPTS) calculates Twiss parameters
   |     and closed orbit coordinates at specified reference points REFPTS.
   | 
   |     Note: REFPTS is an array of increasing indexes that  
   |     select elements from range 1 to length(LATTICE)+1. 
   |     See further explanation of REFPTS in the 'help' for FINDSPOS 
   | 
   |  [TwissData, tune, chrom]  = TWISSRING(...,'chrom', DDP) also calculates
   |     linear dispersion and chromaticity. Dispersion is returned as one 
   |     of the fields in TwissData.
   |     !!! Last argument DDP is a momentum deviation on top
   |     of DP (the second argument) used to calculate and normalize
   |     dispersion and chromaticity. If not supplied
   |     the default value of 1e-8 is used.
   | 
   |     Note: To resolve the integer part of the tune
   |     and the uncertainty of acos(trace(M)/2) it is necessary to
   |     supply sufficient number of REFPTS properly spaced in betatron phase.
   | 
   |  TwisData is a 1-by-REFPTS (1-by-1) structure array with fields
   |        (Some are the same as in the output of LINOPT)
   |        ElemIndex   - integer (element number) in the RING 
   |        SPos        - longitudinal position [m]
   |        ClosedOrbit - closed orbit column vector with 
   |                      components x, px, y, py (momentums, NOT angles)						
   |        Dispersion  - dispersion orbit position 4-by-1 vector with 
   |                      components [eta_x, eta_prime_x, eta_y, eta_prime_y]'
   |                      calculated with respect to the closed orbit with 
   |                      momentum deviation DP
   |        M44         - 4x4 transfer matrix M from the beginning of RING
   |                      to the entrance of the element for specified DP [2]
   |        beta        - [betax, betay] horizontal and vertical Twiss parameter beta
   |        alpha       - [alphax, alphay] horizontal and vertical Twiss parameter alpha
   |        mu          - [mux, muy] horizontal and vertical betatron phase
   |                      !!! NOT 2*PI normalized
   |  
   |  Use MATLAB function CAT to get the data from fields of TwissData into MATLAB arrays.
   |      Example: 
   |      >> TD = twissring(THERING,0,1:length(THERING));
   |      >> BETA = cat(1,TD.beta);
   |      >> S = cat(1,TD.SPos);
   |      >> plot(S,BETA(:,1))
   |   
   |  See also TWISSLINE, LINOPT, TUNECHROM.

.. py:function:: RadIntegrals

   | Wigidx is index of insertion devices
   | 
   |  A.Mash'al, Iranian Light Source Facility, 2018-07-30   

.. py:function:: atx

   | 
   | BEAMDATA=ATX(RING) Computes linear optics, equilibrium emittances and damping times
   | 
   |  RING:     Ring description. If RING is 6D, it is used in OHMIENVELOPE and
   |            a 4D copy may be used for linear optics computation.
   |            If RING is 4D, a 6D copy with default options is used for
   |            OHMIENVELOPE
   | 
   | BEAMDATA=ATX(RING,DP,REFPTS)
   | BEAMDATA=ATX(RING,REFPTS)
   |    Specify the points of interest (Default: 1:length(RING)+1)
   | 
   | BEAMDATA=ATX(RING,DP,...) 
   | BEAMDATA=ATX(RING,...,'dp',DPP)
   |    Specify the momentum deviation (Default: 0)
   | 
   | BEAMDATA=ATX(RING,...,'dct',DCT)
   |    Specify the path lengthening
   | 
   | BEAMDATA=ATX(RING,...,'df',DF)
   |    Specify the RF frequency deviation from nominal
   | 
   | BEAMDATA=ATX(RING,...,'method',OPTICSFUN)
   |    Specify the method for linear optics. Default: @atlinopt6
   |    Allowed values are @atlinopt2, @atlinopt4, @atlinopt6
   | 
   | BEAMDATA=ATX(RING,...,'6d')
   |    By default, linear optics is computed with the 4d version of the lattice.
   |    If method is @atlinopt6 (the default), when specifying '6d' the optics
   |    is computed from the 6d version of the lattice. 
   | 
   | ELEMDATA is a MATLAB structure array as long as the numver of refpoints
   | with fields:
   | 
   |  From atlinopt:
   | 
   |    ElemIndex   - ordinal position in the RING
   |    SPos        - longitudinal position [m]
   |    ClosedOrbit - closed orbit column vector with
   |                  components x, px, y, py (momentums, NOT angles)
   |    Dispersion  - dispersion orbit position vector with
   |                  components eta_x, eta_prime_x, eta_y, eta_prime_y
   |                  calculated with respect to the closed orbit with
   |                  momentum deviation DP
   |    M           - 4x4 transfer matrix M from the beginning of RING
   |                  to the entrance of the element for specified DP [2]
   |    mu          - [ mux, muy] horizontal and vertical betatron phase
   |    beta        - [betax, betay] vector
   |    alpha       - [alphax, alphay] vector
   | 
   |    Other fields depend on the selected optics method, @atlinopt4 or
   |    @atlinopt6
   | 
   |  From ohmienvelope:
   | 
   |    beam66      - 6x6 equilibrium beam matrix
   |    emit66      - 6x6 emittance projections on x-x', y-y' and energy spread
   |    beam44      - intersection of beam66 with dpp=0 (monochromatic beam)
   |    emit44      - emittances of the projections of beam44 on x and y
   |    modemit     - emittance of the 3 eigenmodes
   | 
   | [ELEMDATA,RINGDATA]=ATX(...)  Returns also a structure RINGDATA
   | with fields:
   | 
   |    ll          - Circumference
   |    alpha       - momentum compaction factor
   |    fractunes
   |    fulltunes
   |    nuh         - Tunes
   |    nuv
   |    chromaticity
   |    dampingtime
   |    espread     - Energy spread
   |    blength     - Bunch length
   |    energy
   |    fs          - synchrotron frequency
   |    eloss       - energy loss/turn
   |    synchrophase- synchronous phase
   |    modemittance- Eigen emittances
   |    momcompact  - momentum compaction factor
   | 
   | 
   | The following options are kept for backwards compatibility but are
   | deprecated:
   | 
   | BEAMDATA=ATX(RING,DPP,REFPTS,RADRING,RADINDEX,CAVINDEX)
   |  Radiation must be turned on for emittance computation. This is done by
   |  default using the ATRADON function with default arguments. If this is not
   |  desired, this syntax allows to explicitly enter the radiative lattice
   | 
   | BEAMDATA=ATX(RING,DPP,REFPTS,RADFUNCTION)
   |  RADFUNCTION is substituted to ATRADON to provide the radiative lattice
   |  and indices, in the form:
   |         [RADRING,RADINDEX,CAVINDEX]=RADFUNCTION(RING)
   | 
   |  See also atlinopt atradon ohmienvelope ringpara atsummary

.. py:function:: twissline

   |  
   |  TwissData  = TWISSLINE(LATTICE,DP,TWISSDATAIN) propagates twiss
   |     parameters and closed orbit coordinates from the LINE entrance
   |     given by TWISSDATAIN assuming constant energy deviation DP.
   |     TWISSDATAIN is a 1-by-1 structure with the same field names
   |     as the return argument. (See below)
   |     !!! IMPORTANT: Since  TWISSLINE does not search for closed orbit
   |     its value at the entrance must be supplied in the 
   |     ClosedOrbit field of TWISSDATAIN structure. 
   | 
   |  TwissData  = TWISSLINE(LATTICE,DP,TWISSDATAIN,REFPTS) calculates Twiss parameters
   |     and closed orbit coordinates at specified reference points REFPTS
   | 
   |     Note: REFPTS is an array of increasing indexes that  
   |     select elements from range 1 to length(LATTICE)+1. 
   |     See further explanation of REFPTS in the 'help' for FINDSPOS 
   | 
   |  TwissData  = TWISSLINE(...,'chrom', DDP) also calculates
   |     linear dispersion. Dispersion is returned as one 
   |     of the fields in TwissData.
   |     !!! Last argument DDP is a momentum deviation on top
   |     of DP (the second argument) used to calculate and normalize
   |     dispersion. If not supplied
   |     the default value of 1e-8 is used.
   | 
   |  TwisData is a 1-by-REFPTS (1-by-1 if no REFPTS specified) structure array with fields:
   |        ElemIndex   - integer (element number) in the LINE 
   |        SPos        - longitudinal position [m]
   |        ClosedOrbit - closed orbit column vector with 
   |                      components x, px, y, py (momentums, NOT angles)						
   |        Dispersion  - dispersion orbit position 4-by-1 vector with 
   |                      components [eta_x, eta_prime_x, eta_y, eta_prime_y]'
   |                      calculated with respect to the closed orbit with 
   |                      momentum deviation DP
   |        M44         - 4x4 transfer matrix M from the beginning of LINE
   |                      to the entrance of the element for specified DP [2]
   |        beta        - [betax, betay] horizontal and vertical Twiss parameter beta
   |        alpha       - [alphax, alphay] horizontal and vertical Twiss parameter alpha
   |        mu          - [mux, muy] horizontal and vertical betatron phase
   |                      !!! NOT 2*PI normalized
   |  
   |  Use CAT to get the data from fields of TwissData into MATLAB arrays.
   |      Example: 
   |      >> TD = twissring(THERING,0,1:length(THERING));
   |      >> BETA = cat(1,TD.beta);
   |      >> S = cat(1,TD.SPos);
   |      >> plot(S,BETA(:,1))
   |   
   |  See also TWISSRING, LINOPT, TUNECHROM.

.. py:function:: ringpara

   | 
   | RINGPARA()                 Print parameters of the global variable THERING
   | RINGPARA(RING)             Print parameters of RING. RING may be 4D or 6D
   | 
   | RINGPARA(RING,U0,...)      Supply the total radiation loss in MeV
   | 
   | RINGPARA(...,'dp',DP)      Specify off-momentum For off-momentum lattices,
   |    the equilibrium emittance and energy spread will not be computed
   |    because the quadrupole contribution is missing.
   | 
   | RINGPARA(...,'dct',DCT)    Specify the path lengthening
   | 
   | RINGPARA(...,'df',DF)      Specify the RF frequency deviation
   | 
   | RINGPARAMS=RINGPARA(...)   Return the results in a structure instead of
   |    printing them
   | 
   |   See also ATX ATSUMMARY

