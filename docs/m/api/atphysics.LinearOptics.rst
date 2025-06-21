.. _linearoptics_module:

LinearOptics
============

.. rubric:: Functions


.. list-table::

   * - :func:`amat`
     - find A matrix from one turn map matrix T such that:
   * - :func:`atdampingrates`
     - find tunes and damping rates from one map matrix with radiation
   * - :func:`atlinopt`
     - Performs 4D linear analysis of the COUPLED lattices
   * - :func:`atlinopt2`
     - Performs the linear analysis of UNCOUPLED lattices
   * - :func:`atlinopt4`
     - Performs the 4D linear analysis of COUPLED lattices
   * - :func:`atlinopt6`
     - Performs linear analysis of the lattice
   * - :func:`beam22`
     - computes the beam matrix from the 1-turn transfer matrix
   * - :func:`beam44`
     - computes the coupled beam matrices
   * - :func:`find_betaoids`
     - [H1 H2 H3]=find_betaoids(A)
   * - :func:`find_etaoids`
     - Given the normalizing matrix A, we compute the etaoids
   * - :func:`find_inv_G`
     - This function computes the invariants of a one turn map matrix.
   * - :func:`find_inv_G_fromA`
     - This function computes the invariant matrices G1,G2,G3
   * - :func:`findelemm44`
     - numerically finds the 4x4 transfer matrix of an element
   * - :func:`findelemm66`
     - numerically finds the 6x6 transfer matrix of an element
   * - :func:`findm44`
     - numerically finds the 4x4 transfer matrix of an accelerator lattice
   * - :func:`findm66`
     - numerically finds the 6x6 transfer matrix of an accelerator lattice
   * - :func:`get_dispersion_from_etaoids`
     - get_dispersion_from_etaoids computes dispersion functions (x,px,y,py) at refpts
   * - :func:`jmat`
     - Compute antisymmetric Matrix [O 1; -1 0]
   * - :func:`linopt`
     - performs linear analysis of the COUPLED lattices
   * - :func:`mkSRotationMatrix`
     - (PSI) coordinate transformation matrix
   * - :func:`plotbeta`
     - plots UNCOUPLED! beta-functions
   * - :func:`r_analysis`
     - 

.. py:function:: amat

   | find A matrix from one turn map matrix T such that:
   | 
   |            [Rotx  0    0  ]
   | inv(A)*T*A=[ 0   Rotz  0  ]
   |            [ 0    0   Rots]
   | 
   | Order it so that it is close to the order of x,y,z
   | also ensure that positive modes are before negative so that
   | one has proper symplecticity
   | B. Nash July 18, 2013
   | we order and normalize the vectors via
   |  v_j' jmat(3) v_k = i sgn(j) delta(j,k)

.. py:function:: atdampingrates(m66)

   | find tunes and damping rates from one map matrix with radiation
   | 
   | **[nu,chi]=atdampingrates(m66)**
   | 
   | note that in order to find the damping times, one needs the revolution
   | time T0, then
   | tau1 = T0/chi1, tau2 = T0/chi2, tau3 = T0/chi3

.. py:function:: atlinopt(ring,dp,refpts)

   | Performs 4D linear analysis of the COUPLED lattices
   | 
   | **atlinopt** only works for CONSTANT energy. So no PassMethod should:
   |    1. change the longitudinal momentum dP (cavities, magnets with radiation)
   |    2. have any time dependence (localized impedance, fast kickers etc)
   | 
   | For a more general computation, see ATLINOPT6
   | 
   |  **lindata = atlinopt(ring,dp,refpts)** is a MATLAB structure array with fields
   | 
   |    ElemIndex   - ordinal position in the RING
   |    SPos        - longitudinal position [m]
   |    ClosedOrbit - 4x1 closed orbit vector with components
   |                  x, px, y, py (momentums, NOT angles)
   |    Dispersion  - 4x1 dispersion with components
   |                  eta_x, eta_prime_x, eta_y, eta_prime_y
   |    M44         - 4x4 transfer matrix M from the beginning of RING
   |                  to the entrance of the element [2]
   |    A           - 2x2 matrix A in [3]
   |    B           - 2x2 matrix B in [3]
   |    C           - 2x2 matrix C in [3]
   |    gamma       - gamma parameter of the transformation to eigenmodes
   |    mu          - [mux, muy] horizontal and vertical betatron phase advances
   |    beta        - [betax, betay] vector
   |    alpha       - [alphax, alphay] vector
   | 
   |    All values are specified at the entrance of each element specified in REFPTS.
   |    REFPTS is an array of increasing indexes that  select elements
   |    from the range 1 to length(LINE)+1. Defaults to 1 (initial point)
   |    See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |  **[lindata,nu] = atlinopt()** returns a vector of linear tunes
   |    [nu_u , nu_v] for two normal modes of linear motion [1]
   | 
   |  **[lindata,nu, ksi] = atlinopt() returns a vector of chromaticities ksi = d(nu)/(dp/p)**
   |    [ksi_u , ksi_v] - derivatives of [nu_u , nu_v]
   | 
   |  **[...] = atlinopt(...,'orbit',orbitin)**
   |  **[...] = atlinopt(ring,dp,refpts,orbitin)  (deprecated syntax)**
   |    Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
   |    of initial conditions is used: [x0; px0; y0; py0; DP; 0].
   |    The sixth component is ignored.
   |    This syntax is useful to specify the entrance orbit if RING is not a
   |    ring or to avoid recomputing the closed orbit if is already known.
   | 
   |  **[...] = atlinopt(...,'coupled',flag)**
   |    If flag is false, a faster computation is performed
   |    assuming no coupling in the lattice. Default: true
   | 
   |  **[...] = atlinopt(...,'ct',ct)**
   |    Instead of computing the linear optics for  the specified DP/P,
   |    computes for the path lenghening specified by CT.
   |    The DP argument is ignored.
   | 
   |  **[...] = atlinopt(...,'twiss_in',twissin)**
   |    Computes the optics for a transfer line.
   | 
   |  TWISSIN is a scalar structure with fields:
   |    ClosedOrbit - 4x1 initial closed orbit. Default: zeros(4,1)
   |    Dispersion  - 4x1 initial dispersion.   Default: zeros(4,1)
   |    mu          - [ mux, muy] horizontal and vertical betatron phase
   |    beta        - [betax0, betay0] vector
   |    alpha       - [alphax0, alphay0] vector
   | 
   |  Difference with linopt: Fractional tunes 0<=tune<1
   | 			  Dispersion output
   | 			  Alpha output
   | 			  Phase advance output
   | 			  Option to skip closed orbit search
   |   REFERENCES
   |     [1] D.Edwars,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
   |     [2] E.Courant, H.Snyder
   |     [3] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams, vol.2 (1999)
   | 
   |   See also atlinopt2 atlinopt4 atlinopt6 atx tunechrom

.. py:function:: atlinopt2(ring,refpts)

   | Performs the linear analysis of UNCOUPLED lattices
   | 
   | **[ringdata,elemdata] = atlinopt2(ring,refpts)**
   | 
   |  IMPORTANT!!! **atlinopt2** assumes a constant momentum deviation.
   |    PassMethods used for any element in the RING SHOULD NOT
   |    1.change the longitudinal momentum dP
   |      (cavities , magnets with radiation, ...)
   |    2.have any time dependence (localized impedance, fast kickers, ...)
   | 
   | RINGDATA is a structure array with fields:
   |    tune          1x2 tunes
   |    chromaticity  1x2 chromaticities (only with get_chrom or get_w flags)
   | 
   | ELEMDATA is a structure array with fields:
   |    SPos        - longitudinal position [m]
   |    ClosedOrbit - 4x1 closed orbit vector with components
   |                  x, px, y, py (momentums, NOT angles)
   |    Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
   |    M           - 4x4 transfer matrix M from the beginning of RING
   |                  to the entrance of the element [2]
   |    beta        - [betax, betay] vector
   |    alpha       - [alphax, alphay] vector
   |    mu          - [mux, muy] horizontal and vertical betatron phase advances
   |    W           - [Wx, Wy]  Chromatic amplitude function [3] (only with the
   |                            get_w flag)
   | 
   |    All values are specified at the entrance of each element specified in REFPTS.
   |    REFPTS is an array of increasing indexes that  select elements
   |    from the range 1 to length(LINE)+1. Defaults to 1 (initial point)
   |    See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |    Use the Matlab function "cat" to get the data from fields of ELEMDATA as MATLAB arrays.
   |    Example:
   |    >> **[ringdata, elemdata] = atlinopt2(ring,1:length(ring))**;
   |    >> beta = cat(1,elemdata.beta);
   |    >> s = cat(1,elemdata.SPos);
   |    >> plot(S,beta)
   | 
   |  **[...] = atlinopt2(...,'get_chrom')**
   |    Trigger the computation of chromaticities
   | 
   |  **[...] = atlinopt2(...,'get_w')**
   |    Trigger the computation of chromatic amplitude functions (time consuming)
   | 
   |  **[...] = atlinopt2(...,'dp',dpp)**
   |    Analyses the off-momentum lattice by specifying the central
   |    off-momentum value
   | 
   |  **[...] = atlinopt2(...,'ct',dct)**
   |    Analyses the off-momentum lattice by specifying the path lengthening
   |    corresponding to a frequency shift. The resulting deltap/p is returned
   |    in the 5th component of the ClosedOrbit field.
   | 
   |  **[...] = atlinopt2(...,'orbit',orbitin)**
   |    Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
   |    of initial conditions is used: [x0; px0; y0; py0; DP; 0].
   |    The sixth component is ignored.
   |    This syntax is useful to specify the entrance orbit if RING is not a
   |    ring or to avoid recomputing the closed orbit if is already known.
   | 
   |  **[...] = atlinopt2(...,'twiss_in',twissin)**
   |    Computes the optics for a transfer line.
   | 
   |  TWISSIN is a scalar structure with fields:
   |    ClosedOrbit - 4x1 initial closed orbit. Default: zeros(4,1)
   |    Dispersion  - 4x1 initial dispersion.   Default: zeros(4,1)
   |    mu          - [ mux, muy] horizontal and vertical betatron phase
   |    beta        - [betax0, betay0] vector
   |    alpha       - [alphax0, alphay0] vector
   | 
   |   REFERENCES
   | 	[1] D.Edwards,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
   | 	[2] E.Courant, H.Snyder
   | 	[3] Brian W. Montague Report LEP Note 165, CERN, 1979
   | 
   |   See also atlinopt atlinopt4 atlinopt6 tunechrom

.. py:function:: atlinopt4(ring,refpts)

   | Performs the 4D linear analysis of COUPLED lattices
   | 
   | **[ringdata,elemdata] = atlinopt4(ring,refpts)**
   | 
   |  IMPORTANT!!! **atlinopt4** assumes a constant momentum deviation.
   |    PassMethods used for any element in the RING SHOULD NOT
   |    1.change the longitudinal momentum dP
   |      (cavities , magnets with radiation, ...)
   |    2.have any time dependence (localized impedance, fast kickers, ...)
   | 
   | RINGDATA is a structure array with fields:
   |    tune          1x2 tunes
   |    chromaticity  1x2 chromaticities (only with get_chrom or get_w flags)
   | 
   | ELEMDATA is a structure array with fields:
   |    SPos        - longitudinal position [m]
   |    ClosedOrbit - 4x1 closed orbit vector with components
   |                  x, px, y, py (momentums, NOT angles)
   |    Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
   |    M           - 4x4 transfer matrix M from the beginning of RING
   |                  to the entrance of the element [2]
   |    A           - 2x2 matrix A in [4]
   |    B           - 2x2 matrix B in [4]
   |    C           - 2x2 matrix C in [4]
   |    gamma       - gamma parameter of the transformation to eigenmodes [4]
   |    beta        - [betax, betay] vector
   |    alpha       - [alphax, alphay] vector
   |    mu          - [mux, muy] Betatron phase advances
   |    W           - [Wx, Wy]  Chromatic amplitude function [3] (only with the
   |                            get_w flag)
   | 
   |    Use the Matlab function "cat" to get the data from fields of ELEMDATA as MATLAB arrays.
   |    Example:
   |    >> **[ringdata, elemdata] = atlinopt4(ring,1:length(ring))**;
   |    >> beta = cat(1,elemdata.beta);
   |    >> s = cat(1,elemdata.SPos);
   |    >> plot(S,beta)
   | 
   |    All values are specified at the entrance of each element specified in REFPTS.
   |    REFPTS is an array of increasing indexes that  select elements
   |    from the range 1 to length(LINE)+1. Defaults to 1 (initial point)
   |    See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |  **[...] = atlinopt4(...,'get_chrom')**
   |    Trigger the computation of chromaticities
   | 
   |  **[...] = atlinopt4(...,'get_w')**
   |    Trigger the computation of chromatic amplitude functions (time consuming)
   | 
   |  **[...] = atlinopt4(...,'dp',dpp)**
   |    Analyses the off-momentum lattice by specifying the central
   |    off-momentum value
   | 
   |  **[...] = atlinopt4(...,'ct',dct)**
   |    Analyses the off-momentum lattice by specifying the path lengthening
   |    corresponding to a frequency shift. The resulting deltap/p is returned
   |    in the 5th component of the ClosedOrbit field.
   | 
   |  **[...] = atlinopt4(...,'df',df)**
   |    Analyses the off-momentum lattice by specifying the RF frequency shift.
   |    The resulting deltap/p is returned in the 5th component of the ClosedOrbit field.
   | 
   |  **[...] = atlinopt4(...,'orbit',orbitin)**
   |    Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
   |    of initial conditions is used: [x0; px0; y0; py0; DP; 0].
   |    The sixth component is ignored.
   |    This syntax is useful to specify the entrance orbit if RING is not a
   |    ring or to avoid recomputing the closed orbit if is already known.
   | 
   |  **[...] = atlinopt4(...,'twiss_in',twissin)**
   |    Computes the optics for a transfer line.
   | 
   |  TWISSIN is a scalar structure with fields:
   |    ClosedOrbit - 4x1 initial closed orbit. Default: zeros(4,1)
   |    Dispersion  - 4x1 initial dispersion.   Default: zeros(4,1)
   |    mu          - [ mux, muy] horizontal and vertical betatron phase
   |    beta        - [betax0, betay0] vector
   |    alpha       - [alphax0, alphay0] vector
   | 
   |   REFERENCES
   | 	[1] D.Edwards,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
   | 	[2] E.Courant, H.Snyder
   | 	[3] Brian W. Montague Report LEP Note 165, CERN, 1979
   | 	[4] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams, vol.2 (1999)
   | 
   |   See also atlinopt2 atlinopt6 tunechrom

.. py:function:: atlinopt6(ring,refpts)

   | Performs linear analysis of the lattice
   | 
   |  **[ringdata,elemdata] = atlinopt6(ring,refpts)**
   | 
   | For circular machines, **atlinopt6** analyses
   | the 4x4 1-turn transfer matrix if radiation is OFF, or
   | the 6x6 1-turn transfer matrix if radiation is ON.
   | 
   | For a transfer line, The "twiss_in" intput must contain either:
   |  - a field 'R', as provided by **atlinopt6**, or
   |  - the fields 'beta' and 'alpha', as provided by ATLINOPT and **atlinopt6**
   | 
   | RINGDATA is a structure array with fields:
   |    tune            Fractional tunes
   |    damping_time    Damping times [s]
   |    chromaticity    Chromaticities
   | 
   | ELEMDATA is a structure array with fields:
   |    SPos        - longitudinal position [m]
   |    ClosedOrbit - 6x1 closed orbit vector with components
   |                  x, px, y, py, dp/d, ct (momentums, NOT angles)
   |    Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
   |    R           - DxDx(D/2) R-matrices
   |    A           - DxD A-matrix
   |    M           - DxD transfer matrix M from the beginning of RING
   |    beta        - [betax, betay]                 1x2 beta vector
   |    alpha       - [alphax, alphay]               1x2 alpha vector
   |    mu          - [mux, muy] 	Betatron phase advances
   |    W           - [Wx, Wy]  Chromatic amplitude function [3] (only with the
   |                            get_w flag)
   | 
   |    Use the Matlab function "cat" to get the data from fields of ELEMDATA as MATLAB arrays.
   |    Example:
   |    >> **[ringdata, elemdata] = atlinopt6(ring,1:length(ring))**;
   |    >> beta = cat(1,elemdata.beta);
   |    >> s = cat(1,elemdata.SPos);
   |    >> plot(S,beta)
   | 
   |    All values are specified at the entrance of each element specified in REFPTS.
   |    REFPTS is an array of increasing indexes that  select elements
   |    from the range 1 to length(LINE)+1. Defaults to 1 (initial point)
   |    See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |  **[...] = atlinopt6(...,'get_chrom')**
   |    Trigger the computation of chromaticities
   | 
   |  **[...] = atlinopt6(...,'get_w')**
   |    Trigger the computation of chromatic amplitude functions (time consuming)
   | 
   |  **[...] = atlinopt6(...,'orbit',orbitin)**
   |    Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
   |    of initial conditions is used: [x0; px0; y0; py0; DP; 0].
   |    The sixth component is ignored.
   |    This syntax is useful to specify the entrance orbit if RING is not a
   |    ring or to avoid recomputing the closed orbit if is already known.
   | 
   |  **[...] = atlinopt6(...,'twiss_in',twissin)**
   |    Computes the optics for a transfer line.
   |        TWISSIN is a scalar structure with fields:
   |            R       4x4x2 R-matrices
   |                    or
   |            beta    [betax0, betay0] vector
   |            alpha	[alphax0, alphay0] vector
   | 
   |  **[...] = atlinopt6(...,'dp',dp)**
   |    Specify off-momentum
   | 
   |  **[...] = atlinopt6(...,'dct',dct)**
   |    Specify the path lengthening
   | 
   |  **[...] = atlinopt6(...,'df',df)**
   |    Specify the RF frequency deviation
   | 
   |   REFERENCES
   |    [1] Etienne Forest, Phys. Rev. E 58, 2481 – Published 1 August 1998
   |    [2] Andrzej Wolski, Phys. Rev. ST Accel. Beams 9, 024001 –
   |        Published 3 February 2006
   |    [3] Brian W. Montague Report LEP Note 165, CERN, 1979
   | 
   |   See also atlinopt2 atlinopt4 tunechrom

.. py:function:: beam22(t)

   | computes the beam matrix from the 1-turn transfer matrix
   | %
   | **beam=beam22(t)**
   |  T:    1-turn transfer matrix
   |  BEAM: envelope matrix
   | 
   | **[beam,tune]=beam22(t)**
   |  also returns the tune

.. py:function:: beam44(a,b,c,gamma)

   | computes the coupled beam matrices
   | 
   | **[beam1,beam2]=beam44(a,b,c,gamma)**
   |  A,B,C,gamma: Coupling parameters, see [1]
   |  BEAM1,BEAM2: Eigen modes
   | 
   | **[beam1,beam2]=beam44(lindata)**
   |  LINDATA: structure with fields A,B,C,gamma
   | 
   | **[beam1,beam2,tune1,tune1]=beam44(...)**
   |  also returns the tunes
   | 
   | [1] Sagan, Rubin, "Linear Analysis of Coupled Lattices"
   |     Phys.Rev.Spec.Top. - Accelerators and Beams, vol2, 1999

.. py:function:: find_betaoids

   | [H1 H2 H3]=find_betaoids(A)
   | Given the normalizing matrix A, we compute the betaoids
   | (in Forest's terminology)
   |  these can be related to the invariants

.. py:function:: find_etaoids

   | Given the normalizing matrix A, we compute the etaoids
   | (in terminology of E. Forest)
   |  these can be related to the dispersion

.. py:function:: find_inv_G

   | This function computes the invariants of a one turn map matrix.
   | The resulting invariant matrices G1,G2,G3 satisfy
   |  M^T G_a M = G_a for a=1,2,3
   |  Algorithm from PhD thesis of B. Nash

.. py:function:: find_inv_G_fromA

   | This function computes the invariant matrices G1,G2,G3
   | starting from the A matrix computed via the function amat().
   | The resulting invariant matrices G1,G2,G3 satisfy
   |  M^T G_a M = G_a for a=1,2,3

.. py:function:: findelemm44(elem, methodname)

   | numerically finds the 4x4 transfer matrix of an element
   |   **findelemm44(elem, methodname)**
   |      ELEM          - the element data structure
   |      METHODNAME    - name of the pass-method function
   |                    (default:  ELEM.PassMethod)
   | 
   |   **m66=findelemm44(...,'orbit',orbitin)  (deprecated syntax)**
   |   **m66=findelemm44(elem, methodname, orbitin)**
   |      ORBITIN       - 6x1 phase space coordinates at the entrance
   |                    (default: zeros(6,1))
   |                    The transverse matrix is momentum-dependent,
   |                    the 5-th component of ORBITIN is used as the DP value
   | 
   |   **m66=findelemm44(...,'energy',energy)**
   |      Use ENERGY and ignore the 'Energy' field of elements
   | 
   |   **m66=findelemm44(...,'particle',particle)**
   |      Use PARTICLE (default is relativistic)
   | 
   | 
   |  See also FINDELEMM66

.. py:function:: findelemm66(elem, methodname)

   | numerically finds the 6x6 transfer matrix of an element
   |   **m66=findelemm66(elem, methodname)**
   |      ELEM          - the element data structure
   |      METHODNAME    - name of the pass-method function
   |                    (default:  ELEM.PassMethod)
   | 
   |   **m66=findelemm66(elem, methodname, orbitin)  (deprecated syntax)**
   |   **m66=findelemm66(...,'orbit',orbitin)**
   |      ORBITIN       - 6-by-1 phase space coordinates at the entrance
   |                    (default: zeros(6,1))
   | 
   |   **m66=findelemm66(...,'energy',energy)**
   |      Use ENERGY and ignore the 'Energy' field of elements
   | 
   |   **m66=findelemm66(...,'particle',particle)**
   |      Use PARTICLE (default is relativistic)
   | 
   |  See also FINDELEMM44

.. py:function:: findm44(lattice)

   | numerically finds the 4x4 transfer matrix of an accelerator lattice
   |  for a particle with relative momentum deviation DP
   | 
   |  IMPORTANT!!!
   |  **findm44** gives a wrong result with 6-d rings.
   |  **findm44** assumes constant momentum deviation.
   |    PassMethod used for any element in the LATTICE SHOULD NOT
   |    1.change the longitudinal momentum dP
   |      (cavities , magnets with radiation, ...)
   |    2.have any time dependence (localized impedance, fast kickers, ...)
   | 
   |  **m44 = findm44(lattice)** finds the full one-turn
   |     matrix at the entrance of the first element
   |     !!! With this syntax **findm44** assumes that the LATTICE
   |     is a ring and first finds the closed orbit
   | 
   |  **[m44,t] = findm44(lattice,refpts)** also returns
   |     4x4 transfer matrixes  between entrance of
   |     the first element and each element indexed by REFPTS.
   |     T is 4x4xlength(REFPTS) 3 dimensional array
   |     so that the set of indexes (:,:,i) selects the 4-by-4
   |     matrix at the i-th reference point.
   | 
   |     Note: REFPTS is an array of increasing indexes that
   |     select elements from range 1 to length(LATTICE)+1.
   |     See further explanation of REFPTS in the 'help' for FINDSPOS
   |     When REFPTS= [ 1 2 .. ] the fist point is the entrance of the
   |     first element and T(:,:,1) - identity matrix
   | 
   |     Note: REFPTS is allowed to go 1 point beyond the
   |     number of elements. In this case the last point is
   |     the EXIT of the last element. If LATTICE is a RING
   |     it is also the entrance of the first element
   |     after 1 turn: T(:,:,end) = M
   | 
   |  **[...] = findm44(ring,...,'dp',dp)**
   |  **[...] = findm44(ring,dp,refpts,...)  (deprecated syntax)**
   |    Computes for the off-momentum DP
   | 
   |  **[...] = findm44(ring,...,'dct',dct)**
   |    Computes for the path lenghening specified by CT.
   | 
   |  **[...] = findm44(ring,...,'df',df)**
   |    Computes for a deviation of RF frequency DF
   | 
   |  **[...] = findm44(ring,...,'orbit',orbitin)**
   |  **[...] = findm44(ring,dp,refpts,orbitin)  (deprecated syntax)**
   |    Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
   |    of initial conditions is used: [x0; px0; y0; py0; DP; 0].
   |    The sixth component is ignored.
   |    This syntax is useful to specify the entrance orbit if RING is not a
   |    ring or to avoid recomputing the closed orbit if is already known.
   | 
   |  **[...] = findm44(...,'full')**
   |    Same as above except that matrices returned in T are full 1-turn
   |    matrices at the entrance of each element indexed by REFPTS.
   | 
   |  **[m44,t,orbit] = findm44(...)**
   |    In addition returns the closed orbit at the entrance of each element
   |    indexed by REFPTS.
   | 
   |  See also FINDM66, FINDORBIT4, ATENABLE_6D, ATDISABLE_6D, CHECK_6D

.. py:function:: findm66(ring)

   | numerically finds the 6x6 transfer matrix of an accelerator lattice
   |   by differentiation of LINEPASS near the closed orbit
   |   **findm66** uses FINDORBIT6 to search for the closed orbit in 6-d
   |   In order for this to work the ring MUST have a CAVITY element
   | 
   |  **m66 = findm66(ring)** finds the full one-turn 6-by-6
   |     matrix at the entrance of the first element
   | 
   | **[...]=findm66(ring,...,'dp',dp)** Specify the momentum deviation when
   |    radiation is OFF (default: 0)
   | 
   | **[...]=findm66(ring,...,'dct',dct)** Specify the path lengthening when
   |    radiation is OFF (default: 0)
   | 
   | **[...]=findm66(ring,...,'df',df)** Specify the RF frequency deviation
   |    radiation is OFF (default: 0)
   | 
   | **[...]=findm66(ring,...,'orbit',orbit)** Specify the orbit at the entrance
   |    of the ring, if known.
   | 
   |  **[m66,t] = findm66(ring,refpts)** in addition to M finds
   |     6-by-6 transfer matrixes  between entrances of
   |     the first element and each element indexed by REFPTS.
   |     T is 6-by-6-by-length(REFPTS) 3 dimentional array.
   | 
   |     REFPTS is an array of increasing indexes that  select elements
   |     from the range 1 to length(RING)+1.
   |     See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |     Note:
   |     When REFPTS= [ 1 2 .. ] the first point is the entrance of the first element
   |     and T(:,:,1) - identity matrix
   |     When REFPTS= [  .. length(RING)+1] the last point is the exit of the last element
   |     and the entrance of the first element after 1 turn: T(:,:, ) = M
   | 
   |  **[...] = findm66(ring,refpts,orbitin)    (deprecated syntax)**
   |  **[...] = findm66(...,'orbit',orbitin)**
   |    Do not search for closed orbit. This syntax is useful to avoid
   |    recomputing the closed orbit if it is already known.
   | 
   |  **[m66,t,orbit] = findm66(...)**
   |    In addition returns the closed orbit at the entrance of each element
   |    indexed by REFPTS.
   | 
   | 
   |  See also FINDM44, FINDORBIT6

.. py:function:: get_dispersion_from_etaoids

   | get_dispersion_from_etaoids computes dispersion functions (x,px,y,py) at refpts
   | using the etaoids (E. Forest terminology) which are computed from the one turn map

.. py:function:: jmat

   | Compute antisymmetric Matrix [O 1; -1 0]
   | 
   |   INPUTS
   |     1. dim - 1,2,3 Dimension of the sqare matrix
   | 
   |   OUPUTS
   |     2. mat - Antisymmetric block Matrix [O 1; -1 0]
   | 
   |   See also symplectify

.. py:function:: linopt(ring,dp,refpts)

   | performs linear analysis of the COUPLED lattices
   |    Notation is the same as in reference [3]
   | 
   | 
   |  **lindata = linopt(ring,dp,refpts)** is a MATLAB structure array with fields
   | 
   |    ElemIndex   - ordinal position in the RING
   |    SPos        - longitudinal position [m]
   |    ClosedOrbit - closed orbit column vector with
   |                  components x, px, y, py (momentums, NOT angles)
   |    Dispersion  - dispersion orbit position vector with
   |                  components eta_x, eta_prime_x, eta_y, eta_prime_y
   |                  calculated with respect to the closed orbit with
   |                  momentum deviation DP
   |    M44         - 4x4 transfer matrix M from the beginning of RING
   |                  to the entrance of the element for specified DP [2]
   |    A           - 2x2 matrix A in [3]
   |    B           - 2x2 matrix B in [3]
   |    C           - 2x2 matrix C in [3]
   |    gamma       - gamma parameter of the transformation to eigenmodes
   |    mu          - [ mux, muy] horizontal and vertical betatron phase
   |    beta        - [betax, betay] vector
   | 
   |    All values are specified at the entrance of each element specified in REFPTS.
   |    REFPTS is an array of increasing indexes that  select elements
   |    from the range 1 to length(LINE)+1.
   |    See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |  **[lindata,nu] = linopt()** returns a vector of linear tunes
   |    [nu_u , nu_v] for two normal modes of linear motion [1]
   | 
   |  **[lindata,nu, ksi] = linopt() returns a vector of chromaticities ksi = d(nu)/(dp/p)**
   |    [ksi_u , ksi_v] - derivatives of [nu_u , nu_v]
   | 
   |  See also FINDSPOS TWISSRING TUNECHROM
   | 
   |    [1] D.Edwars,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
   |    [2] E.Courant, H.Snyder
   |    [3] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams, vol.2 (1999)

.. py:function:: mkSRotationMatrix

   | (PSI) coordinate transformation matrix
   |  that describes s-rotation of the ELEMENT
   | 
   |  |   cos(psi)     0       sin(psi)     0         0       0     |
   |  |      0     cos(psi)        0      sin(psi)    0       0     |
   |  |  -sin(psi)     0       cos(psi)     0         0       0     |
   |  |      0     -sin(psi)       0      cos(psi)    0       0     |
   |  |      0         0           0        0         1       0     |
   |  |      0         0           0        0         0       1     |
   | 
   |  Note: AT defines 'positive' s-rotation as clockwise,
   |        looking in the dirction of the beamm
   | 

.. py:function:: plotbeta(ring)

   | plots UNCOUPLED! beta-functions
   |  **plotbeta(ring)** calculates beta functions of the lattice RING
   |  **plotbeta** with no argumnts uses THERING as the default lattice
   |   Note: **plotbeta** uses FINDORBIT4 and LINOPT which assume a lattice
   |   with NO accelerating cavities and NO radiation
   | 
   |  See also PLOTCOD

.. py:function:: r_analysis


