Orbit
=====

.. rubric:: Functions


xorbit_ct XORBIT_CT	Private function used by findsyncorbit
findorbit4 FINDORBIT4 finds closed orbit in the 4-d transverse phase
plotcod PLOTCOD Closed Orbit Distortion
xorbit_6 XORBIT_6	Private function used by findorbit6
findsyncorbit FINDSYNCORBIT finds closed orbit, synchronous with the RF cavity
xorbit_dp XORBIT_DP  Private function used by findorbit4
findorbit FINDORBIT find the closed orbit
findorbit6 FINDORBIT6 finds closed orbit in the full 6-d phase space

.. py:function:: xorbit_ct

   | ORBIT=XORBIT_CT(RING,DCT)	Closed orbit for fixed path lengthening.

.. py:function:: findorbit4

   |  space by numerically solving  for a fixed point of the one turn
   |  map M calculated with LINEPASS
   | 
   |          (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)
   | 
   |     under the CONSTANT MOMENTUM constraint, dP2 = dP1 = dP and
   |     there is NO constraint on the 6-th coordinate CT
   | 
   |  IMPORTANT !!!
   |  FINDORBIT4 gives a wrong result with 6-d rings. If you still want to
   |     get the result, set 'strict' to -1 or 0.
   |  FINDORBIT4 imposes a constraint on dP and relaxes
   |     the constraint on the revolution frequency. A physical storage
   |     ring does exactly the opposite: the momentum deviation of a
   |     particle on the closed orbit settles at the value
   |     such that the revolution is synchronous with the RF cavity
   | 
   |                  HarmNumber*Frev = Frf
   | 
   |     To impose this artificial constraint in FINDORBIT4
   |     PassMethod used for any elemen SHOULD NOT
   |     1. change the longitudinal momentum dP (cavities , magnets with radiation)
   |     2. have any time dependence (localized impedance, fast kickers etc)
   | 
   | FINDORBIT4(RING) is 4x1 vector - fixed point at the
   |     entrance of the 1-st element of the RING (x,px,y,py)
   | 
   | FINDORBIT4(RING,REFPTS) is 4xLength(REFPTS)
   |    array of column vectors - fixed points (x,px,y,py)
   |    at the entrance of each element indexed by the REFPTS array.
   |    REFPTS is an array of increasing indexes that  select elements
   |    from the range 1 to length(RING)+1.
   |    See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   | FINDORBIT4(RING,DP,REFPTS,...) Obsolete syntax
   | FINDORBIT4(RING, ...,'strict',STRICT) Default: 0.
   |    If STRICT is -1, check_6d is skipped
   |    If STRICT is  0, check_6d emits a warning for non-6d rings.
   |    If STRICT is  1, check_6d emits an error for non-6d rings.
   | 
   | FINDORBIT4(RING,...,'dp',DP)   Specifies the off-momentum
   | 
   | FINDORBIT4(RING,...,'dct',DCT) Specifies the path lengthening
   | 
   | FINDORBIT4(RING,...,'df',DF) Specifies RF frequency deviation
   | 
   | FINDORBIT4(RING,dP,REFPTS,GUESS)
   | FINDORBIT4(...,'guess',GUESS)     The search for the fixed point
   |    starts from initial condition GUESS. Otherwise the search starts from
   |    [0; 0; 0; 0; 0; 0]. GUESS must be a 6x1 vector.
   | 
   | FINDORBIT4(...,'orbit',ORBIT)     Specify the orbit at the entrance
   |    of the ring, if known. FINDORBIT4 will then transfer it to the
   |    reference points. ORBIT must be a 6x1 vector.
   | 
   | [ORBIT, FIXEDPOINT] = FINDORBIT4( ... )
   | 	The optional second return parameter is a 6x1 vector:
   |    closed orbit at the entrance of the RING.
   | 
   |  See also FINDSYNCORBIT, FINDORBIT6, ATDISABLE_6D, CHECK_6D

.. py:function:: plotcod

   |  PLOTCOD(RING,DP) finds and plots closed orbit for a given momentum 
   |   deviation DP. It calls FINDORBIT4 which assumes a lattice
   |   with NO accelerating cavities and NO radiation

.. py:function:: xorbit_6

   | ORBIT=XORBIT_6(RING)	Closed orbit for 6D motion

.. py:function:: findsyncorbit

   |  and momentum deviation dP (first 5 components of the phase space vector)
   |  by numerically solving  for a fixed point
   |  of the one turn map M calculated with LINEPASS
   | 
   |        (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)
   | 
   |     under constraints CT2 - CT1 =  dCT = C(1/Frev - 1/Frev0) and dP2 = dP1 , where
   |     Frev0 = Frf0/HarmNumber is the design revolution frequency
   |     Frev  = (Frf0 + dFrf)/HarmNumber is the imposed revolution frequency
   | 
   |  IMPORTANT!!!
   |  FINDSYNCORBIT gives a wrong result with 6-d rings. If you still want to
   |     get the result, set 'strict' to -1 or 0.
   |  FINDSYNCORBIT imposes a constraint (CT2 - CT1) and
   |     dP2 = dP1 but no constraint on the value of dP1, dP2
   |     The algorithm assumes time-independent fixed-momentum ring
   |     to reduce the dimensionality of the problem.
   |     To impose this artificial constraint in FINDSYNCORBIT
   |     PassMethod used for any element SHOULD NOT
   |     1. change the longitudinal momentum dP (cavities , magnets with radiation)
   |     2. have any time dependence (localized impedance, fast kickers etc).
   | 
   | 
   |  FINDSYNCORBIT(RING) is 5x1 vector - fixed point at the
   | 		entrance of the 1-st element of the RING (x,px,y,py,dP)
   | 
   |  FINDSYNCORBIT(RING,REFPTS) is 5xLength(REFPTS)
   |    array of column vectors - fixed points (x,px,y,py,dP)
   |    at the entrance of each element indexed by the REFPTS array.
   |    REFPTS is an array of increasing indexes that  select elements
   |    from the range 1 to length(RING)+1.
   |    See further explanation of REFPTS in the 'help' for FINDSPOS
   |  FINDSYNCORBIT(RING, ...,'strict',STRICT) Default: 0.
   |    If STRICT is -1, check_6d is skipped
   |    If STRICT is  0, check_6d emits a warning for non-6d rings.
   |    If STRICT is  1, check_6d emits an error for non-6d rings.
   | 
   |  See also FINDSYNCORBIT, FINDORBIT6, ATDISABLE_6D, CHECK_6D

.. py:function:: xorbit_dp

   | ORBIT=XORBIT_DP(RING,DP)   Closed orbit for fixed off-momentum.

.. py:function:: findorbit

   | 
   | Depending on the lattice, FINDORBIT will:
   |  - use findorbit6 if radiation is ON,
   |  - use findsyncorbit if radiation is OFF and ct is specified,
   |  - use findorbit4 otherwise
   | 
   | [ORBIT,FIXEDPOINT]=FINDORBIT(RING,REFPTS)
   | 
   |  ORBIT:        6xNrefs, closed orbit at the selected reference points
   |  FIXEDPOINT:   6x1 closed orbit at the entrance of the lattice
   | 
   | [...]=FINDORBIT(RING,...,'dp',DP) Specify the momentum deviation when
   |    radiation is OFF (default: 0)
   | 
   | [...]=FINDORBIT(RING,...,'dct',DCT) Specify the path lengthening when
   |    radiation is OFF (default: 0)
   | 
   | [...]=FINDORBIT(RING,...,'df',DF) Specify the RF frequency deviation
   |    radiation is OFF (default: 0)
   | 
   | [...]=FINDORBIT(RING,...,'orbit',ORBIT) Specify the orbit at the entrance
   |    of the ring, if known. FINDORBIT will then transfer it to the reference points
   | 
   | [...]=FINDORBIT(RING,...,'guess',GUESS) The search will start with GUESS as
   |    initial conditions
   | 
   | For other keywords, see the underlying functions.
   | 
   |  See also FINDORBIT4, FINDSYNCORBIT, FINDORBIT6.

.. py:function:: findorbit6

   |  by numerically solving  for a fixed point of the one turn
   |  map M calculated with RINGPASS
   | 
   |  (X, PX, Y, PY, DP, CT2 ) = M (X, PX, Y, PY, DP, CT1)
   | 
   |  with constraint % CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)
   | 
   |  IMPORTANT!!! FINDORBIT6 is a realistic simulation
   |  1. The Frf frequency in the RF cavities (may be different from Frf0)
   |     imposes the synchronous condition
   |     CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)
   |  2. The algorithm numerically calculates
   |     6-by-6 Jacobian matrix J6. In order for (J-E) matrix
   |     to be non-singular it is NECESSARY to use a realistic
   |     PassMethod for cavities with non-zero momentum kick
   |     (such as RFCavityPass).
   |  3. FINDORBIT6 can find orbits with radiation.
   |     In order for the solution to exist the cavity must supply
   |     adequate energy compensation.
   |     In the simplest case of a single cavity, it must have
   |     'Voltage' field set so that Voltage > Erad - energy loss per turn
   |  4. FINDORBIT6 starts the search from [ 0 0 0 0 0 0 ]', unless
   |     the third argument is specified: FINDORBIT6(RING,REFPTS,GUESS)
   |     There exist a family of solutions that correspond to different RF buckets
   |     They differ in the 6-th coordinate by C*Nb/Frf. Nb = 1 .. HarmNum-1
   |  5. The value of the 6-th coordinate found at the cavity gives
   |     the equilibrium RF phase. If there is no radiation the phase is 0;
   | 
   |  FINDORBIT6(RING) is 6x1 vector - fixed point at the
   | 		entrance of the 1-st element of the RING (x,px,y,py,dp,ct)
   | 
   |  FINDORBIT6(RING,REFPTS) is 6xLength(REFPTS)
   |    array of column vectors - fixed points (x,px,y,py,dp,ct)
   |    at the entrance of each element indexed by the REFPTS array.
   |    REFPTS is an array of increasing indexes that  select elements
   |    from the range 1 to length(RING)+1.
   |    See further explanation of REFPTS in the 'help' for FINDSPOS
   | 
   |  FINDORBIT6(...,'dp',DP)
   |    Specify the off-momentum. The RF frequency will be adjusted to get the
   |    desired value
   | 
   |  FINDORBIT6(...,'dct',DCT)
   |    Specify the path lengthening. The RF frequency will be adjusted to get
   |    the desired value
   | 
   |  FINDORBIT6(...,'df',DF)
   |    Specify the RF frequency deviation
   | 
   |  FINDORBIT6(RING,REFPTS,GUESS)
   |  FINDORBIT6(...,'guess',GUESS)     The search for the fixed point
   | 	starts from initial condition GUESS. Otherwise the search starts from
   |    the synchronous phase. GUESS must be a 6x1 vector.
   | 
   |  FINDORBIT6(...,'orbit',ORBIT)     Specify the orbit at the entrance
   |    of the ring, if known. FINDORBIT6 will then transfer it to the
   |    reference points. ORBIT must be a 6x1 vector.
   | 
   |  [ORBIT, FIXEDPOINT] = FINDORBIT6( ... )
   | 	The optional second return parameter is a 6x1 vector:
   |    closed orbit at the entrance of the RING.
   | 
   |  See also FINDORBIT4, FINDSYNCORBIT.

