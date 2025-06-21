.. _parametersummaryfunctions_module:

ParameterSummaryFunctions
=========================

.. rubric:: Functions


.. list-table::

   * - :func:`RadIntegrals`
     - calcuate the contribution to the radiation integrals of a Wiggler.
   * - :func:`atsummary`
     - Print out various parameters of the current AT lattice
   * - :func:`atx`
     - Computes and displays global information
   * - :func:`ringpara`
     - Print out various parameters of the current AT lattice
   * - :func:`twissline`
     - calculates linear optics functions for an UNCOUPLED transport line
   * - :func:`twissring`
     - calculates linear optics functions for an UNCOUPLED ring

.. py:function:: RadIntegrals

   | calcuate the contribution to the radiation integrals of a Wiggler.
   | Wigidx is index of insertion devices
   | 
   |  A.Mash'al, Iranian Light Source Facility, 2018-07-30

.. py:function:: atsummary()

   | Print out various parameters of the current AT lattice
   | 
   | **atsummary()**                Print parameters of the global variable THERING
   | **atsummary(ring)**            Print parameters of RING. RING may be 4D or 6D
   | 
   | **atsummary(...,'dp',dp)**     Specify off-momentum. For off-momentum lattices,
   |    the equilibrium emittance and energy spread will not be computed
   |    because the quadrupole contribution is missing.
   | 
   | **atsummary(...,'dct',dct)**   Specify the path lengthening
   | 
   | **atsummary(...,'df',df)**     Specify the RF frequency deviation
   | 
   | **ringparams=atsummary(...)**  Return the results in a structure instead of
   |    printing them
   | 
   | **ringparams=atsummary(...,'display')**  Restore the printout of parameters
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

.. py:function:: atx(ring)

   | Computes and displays global information
   | 
   | **beamdata=atx(ring)** Computes linear optics, equilibrium emittances and damping times
   | 
   |  RING:     Ring description. If RING is 6D, it is used in OHMIENVELOPE and
   |            a 4D copy may be used for linear optics computation.
   |            If RING is 4D, a 6D copy with default options is used for
   |            OHMIENVELOPE
   | 
   | **beamdata=atx(ring,dp,refpts)**
   | **beamdata=atx(ring,refpts)**
   |    Specify the points of interest (Default: 1:length(RING)+1)
   | 
   | **beamdata=atx(ring,dp,...)**
   | **beamdata=atx(ring,...,'dp',dpp)**
   |    Specify the momentum deviation (Default: 0)
   | 
   | **beamdata=atx(ring,...,'dct',dct)**
   |    Specify the path lengthening
   | 
   | **beamdata=atx(ring,...,'df',df)**
   |    Specify the RF frequency deviation from nominal
   | 
   | **beamdata=atx(ring,...,'method',opticsfun)**
   |    Specify the method for linear optics. Default: @atlinopt6
   |    Allowed values are @atlinopt2, @atlinopt4, @atlinopt6
   | 
   | **beamdata=atx(ring,...,'6d')**
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
   | **[elemdata,ringdata]=atx(...)**  Returns also a structure RINGDATA
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
   | **beamdata=atx(ring,dpp,refpts,radring,radindex,cavindex)**
   |  Radiation must be turned on for emittance computation. This is done by
   |  default using the ATRADON function with default arguments. If this is not
   |  desired, this syntax allows to explicitly enter the radiative lattice
   | 
   | **beamdata=atx(ring,dpp,refpts,radfunction)**
   |  RADFUNCTION is substituted to ATRADON to provide the radiative lattice
   |  and indices, in the form:
   |         [RADRING,RADINDEX,CAVINDEX]=RADFUNCTION(RING)
   | 
   |  See also atlinopt atradon ohmienvelope ringpara atsummary

.. py:function:: ringpara()

   | Print out various parameters of the current AT lattice
   | 
   | **ringpara()**                 Print parameters of the global variable THERING
   | **ringpara(ring)**             Print parameters of RING. RING may be 4D or 6D
   | 
   | **ringpara(ring,u0,...)**      Supply the total radiation loss in MeV
   | 
   | **ringpara(...,'dp',dp)**      Specify off-momentum For off-momentum lattices,
   |    the equilibrium emittance and energy spread will not be computed
   |    because the quadrupole contribution is missing.
   | 
   | **ringpara(...,'dct',dct)**    Specify the path lengthening
   | 
   | **ringpara(...,'df',df)**      Specify the RF frequency deviation
   | 
   | **ringparams=ringpara(...)**   Return the results in a structure instead of
   |    printing them
   | 
   |   See also ATX ATSUMMARY

.. py:function:: twissline(lattice,dp,twissdatain)

   | calculates linear optics functions for an UNCOUPLED transport line
   | 
   |  **twissdata  = twissline(lattice,dp,twissdatain)** propagates twiss
   |     parameters and closed orbit coordinates from the LINE entrance
   |     given by TWISSDATAIN assuming constant energy deviation DP.
   |     TWISSDATAIN is a 1-by-1 structure with the same field names
   |     as the return argument. (See below)
   |     !!! IMPORTANT: Since  **twissline** does not search for closed orbit
   |     its value at the entrance must be supplied in the
   |     ClosedOrbit field of TWISSDATAIN structure.
   | 
   |  **twissdata  = twissline(lattice,dp,twissdatain,refpts)** calculates Twiss parameters
   |     and closed orbit coordinates at specified reference points REFPTS
   | 
   |     Note: REFPTS is an array of increasing indexes that
   |     select elements from range 1 to length(LATTICE)+1.
   |     See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |  **twissdata  = twissline(...,'chrom', ddp)** also calculates
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

.. py:function:: twissring(lattice,dp)

   | calculates linear optics functions for an UNCOUPLED ring
   | 
   |  **[twissdata, tune]  = twissring(lattice,dp)** calculates twiss parameters
   |     and closed orbit coordinates at the RING entrance assuming
   |     constant energy deviation DP.
   | 
   |  **[twissdata, tune]  = twissring(lattice,dp,refpts)** calculates Twiss parameters
   |     and closed orbit coordinates at specified reference points REFPTS.
   | 
   |     Note: REFPTS is an array of increasing indexes that
   |     select elements from range 1 to length(LATTICE)+1.
   |     See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |  **[twissdata, tune, chrom]  = twissring(...,'chrom', ddp)** also calculates
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
   |      >> **td = twissring(thering,0,1:length(thering))**;
   |      >> BETA = cat(1,TD.beta);
   |      >> S = cat(1,TD.SPos);
   |      >> plot(S,BETA(:,1))
   | 
   |  See also TWISSLINE, LINOPT, TUNECHROM.

