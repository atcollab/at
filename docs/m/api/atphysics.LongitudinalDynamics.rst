.. _longitudinaldynamics_module:

LongitudinalDynamics
====================

.. rubric:: Functions


.. list-table::

   * - :func:`cavityoff`
     - CAVITYOFF turns cavities OFF
   * - :func:`phis`
     - phase = phis(U0MeV,VrfMV)
   * - :func:`atRFacc`
     - ATRFACC Computes RF acceptance of the ring
   * - :func:`atSetCavityPhase`
     - SETCAVITYPHASE     Set the TimeLag attribute of RF cavities
   * - :func:`atBunchLength`
     - bunch length due to the potential well effect
   * - :func:`mcf`
     - MCF momentum compaction factor
   * - :func:`cavityon`
     - CAVITYON turns Cavities ON
   * - :func:`nus`
     - NUS Computes synchrotron tune from RF parameters
   * - :func:`RFacc`
     - RFACC Computes the RF acceptance with linear formula
   * - :func:`BunchLength`
     - 
   * - :func:`atsetcavity`
     - ATSECAVITY Set the parameters of RF cavities

.. py:function:: cavityoff

   |   Sets PassMethod to DriftPass or IdentityPass depending
   |   on the value of 'Length' field
   | 
   |   See also CAVITYON, RADIATIONON, RADIATIONOFF, SETCAVITY

.. py:function:: phis

   | 
   |  this function returns the synchronous phase in radians
   |  input:
   |  U0MeV is energy loss per turn in MeV
   |  VrfMV is the RF voltage in MV

.. py:function:: atRFacc

   |  delta_max_rf = atRFacc(ring)
   |    The functions computes the RF acceptance of the ring
   |    ring is tha at lattice without radiation
   |    delta_max_rf is the RF acceptance
   |    
   |   See also RFacc

.. py:function:: atSetCavityPhase

   | 
   | NEWRING=SETCAVITYPHASE(RING)
   |    Adjust the TimeLag attribute of RF cavities based on frequency,
   |    voltage and energy loss per turn, so that the synchronous phase is zero.
   |    An error occurs if all cavities do not have the same frequency.
   | 
   | NEWRING=SETCAVITYPHASE(...,'refpts',CAVPTS)
   |    CAVPTS is the location of RF cavities. This allows to ignore harmonic
   |    cavities.
   | 
   | WARNING: This function modifies the time reference,
   | this should be avoided
   | 
   | NEWRING=SETCAVITYPHASE(...,'method',METHOD)
   |    Choose the method for computing the energy loss per turn
   | 
   |  METHOD:   'integral': (default) The losses are obtained from
   |                        Losses = Cgamma / 2pi * EGeV^4 * I2
   |                        Takes into account bending magnets and wigglers.
   |            'tracking': The losses are obtained by tracking without cavities.
   |                        Needs radiation ON, takes into account all radiating elements.

.. py:function:: atBunchLength

   |  the output is the zerocurrent bunch length x bunch lengthening
   | 
   |    BL = atBunchLength(ring,Ib,Zn)
   | 
   |  Ib is the bunch current [A] (it may be a vector for multiple values)
   |  Zn is the longitudinal broadband impedance [Ohms]
   |  ring is the at ring without radiation
   |  BL is the bunch length in metres 
   | 
   |    see also: BunchLength

.. py:function:: mcf

   |  MCF(RING) calculates momentum compaction factor of RING
   | 
   |  MCF(RING,DPP) computes the momentum compaction for off-momentum DPP
   | 
   |  IMPORTANT!!!
   |  MCF gives a wrong result with 6-d rings. The RING should be set to 4d.
   |  See also: ATDISABLE_6D, CHECK_6D

.. py:function:: cavityon

   | 
   |  CAVITYON looks for elements that have field Frequency
   |     and sets PassMethod for them to RFCavityPass
   |  CAVITYON(ENERGY)
   |     In addition sets the E0 field of the global variable GLOBVAL
   |     to energy - design energy [eV]
   |     If GLOBVAL does not exist CAVITYON creates it
   | 
   |  See also CAVITYOFF, RADIATIONON, RADIATIONOFF, SETCAVITY

.. py:function:: nus

   |  Nus = nus (VrfMV, alpha, U0MeV, E0MeV, h)
   |  this function return the synchrotron tune
   |  input:
   |  VrfMV is the RF voltage in MV
   |  alpha is the momentum compaction factor
   |  U0MeV is the energy lost per turn in MeV
   |  E0MeV is the beam energy in MeV
   |  h is the harmonic number

.. py:function:: RFacc

   |    delta_max_rf = RFacc(Vrf,U0,E0,h,alpha)
   | 
   |    This function computes the RF acceptance
   |    Vrf is the RF voltage in V
   |    U0 is the energy loss per turn in eV
   |    E0 is the energy of the beam in eV
   |    h is the harmonic number
   |    alpha is the momentum compaction factor
   | 
   |   See also atRFacc

.. py:function:: BunchLength

   |  bunch length due to the potential well effect
   |  the output is the zerocurrent bunch length x bunch lengthening
   | 
   |    BL = BunchLength (Ib,Zn,Vrf,U0,E0,h,alpha,sigdelta,circ)
   | 
   |  Ib is the bunch current [A] (it may be a vector for multiple values)
   |  Zn is the longitudinal broadband impedance [Ohms]
   |  Vrf is the RF voltage [V] (it may be a vector for multiple values)
   |  U0 is the energy loss around the ring [eV]
   |  E0 is the beam energy [eV]
   |  h is the harmonic number
   |  alpha is the momentum compaction factor
   |  sigmadelta is the energy spread
   |  circ is the ring circumference
   |  
   |    see also: atBunchLength

.. py:function:: atsetcavity

   | 
   | ATSETCAVITY may be used in two modes:
   | 
   | Upgrade mode
   | ===================================================
   |  By default, ATSETCAVITY will act on the "main" cavities: they are defined by the
   |  cavpts ring property, or if absent by cavities at the lowest frequency.
   | 
   | NEWRING=ATSETCAVITY(RING,...,'Frequency',FREQUENCY,...)
   |    Set the cavity frequency [Hz]. FREQUENCY is a scalar or an array as
   |    long as the list of selected cavities
   | 
   | NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal',...)
   |    Set the cavity frequency to the nominal value according to
   |    circumference and harmonic number
   | 
   | NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal','dp',DP)
   |    Set the cavity frequency to the nominal value for the specified dp
   | 
   | NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal','dct',DCT)
   |    Set the cavity frequency to the nominal value for the specified dct
   | 
   | NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal','df',DF)
   |    Set the cavity frequency to the nominal value + df
   | 
   | NEWRING=ATSETCAVITY(RING,...,'Voltage',VOLTAGE,...)
   |    Set the total voltage (all cells) [V]. VOLTAGE will be distributed over the
   |    cells: CELL_VOLTAGE = VOLTAGE / PERIODICITY.
   |    Then if CELL_VOLTAGE is a scalar, it will be equally shared among the
   |    selected cavities. Otherwise it is an array as long as the list of
   |    selected cavities.
   | 
   | NEWRING=ATSETCAVITY(RING,...,'HarmNumber',H,...)
   |    Set the harmonic number. H is a scalar or an array as
   |    long as the list of selected cavities
   | 
   | NEWRING=ATSETCAVITY(RING,...,'TimeLag',TIMELAG,...)
   |    Set the time lag [m], . TIMELAG is a scalar or an array as
   |    long as the list of selected cavities
   | 
   | NEWRING=ATSETCAVITY(RING,...,'cavpts',CAVPTS)
   |    CAVPTS is the location of the selected RF cavities. The default is to act on the
   |    "main" cavities: they are defined by the cavpts ring property, or if absent by
   |    cavities at the lowest frequency.
   | 
   |   NOTES
   |   1. In this mode, the radiation state of the lattice is not modified.
   | 
   | 
   | Compatibility mode
   | ===================================================
   | NEWRING = ATSETCAVITY(RING,RFV,RADFLAG,HARM_NUMBER)
   |   RING         Ring structure
   |   RFV          RF voltage (full ring) [V]
   |   RADFLAG      0/1: activate/desactivate radiation (atradon/atradoff)
   |   HARMNUMBER 	Harmonic number (full ring)
   | 
   |   NOTES
   |   1. This mode is deprecated and should be replaced by
   |        RING=ATSETCAVITY(RING,'Frequency','nominal','HarmNumber',HARM_NUMBER, 'Voltage',RFV)
   |        RING=atSetCavityPhase(RING) (optional)
   |        RING=atenable_6d(RING)      (optional)
   |   2. All the N cavities will have a voltage RFV/N
   |   3. sets the synchronous phase of the cavity assuming radiation is turned
   |      on radflag says whether or not we want radiation on, which affects
   |      synchronous phase.
   | 
   |   See also atSetCavityPhase, atsetRFcavity, atenable_6d, atdisable_6d, atgetU0

