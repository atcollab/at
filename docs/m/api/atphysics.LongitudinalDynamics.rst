.. _longitudinaldynamics_module:

LongitudinalDynamics
====================

.. rubric:: Functions


.. list-table::

   * - :func:`BunchLength`
     - 
   * - :func:`RFacc`
     - Computes the RF acceptance with linear formula
   * - :func:`atBunchLength`
     - bunch length due to the potential well effect
   * - :func:`atRFacc`
     - Computes RF acceptance of the ring
   * - :func:`atSetCavityPhase`
     - SETCAVITYPHASE     Set the TimeLag attribute of RF cavities
   * - :func:`atsetcavity`
     - ATSECAVITY Set the parameters of RF cavities
   * - :func:`cavityon`
     - turns Cavities ON
   * - :func:`mcf`
     - momentum compaction factor
   * - :func:`nus`
     - Computes synchrotron tune from RF parameters
   * - :func:`phis`
     - phase = phis(U0MeV,VrfMV)

.. py:function:: BunchLength


.. py:function:: RFacc(vrf,u0,e0,h,alpha)

   | Computes the RF acceptance with linear formula
   |    **delta_max_rf = RFacc(vrf,u0,e0,h,alpha)**
   
   |    This function computes the RF acceptance
   |    Vrf is the RF voltage in V
   |    U0 is the energy loss per turn in eV
   |    E0 is the energy of the beam in eV
   |    h is the harmonic number
   |    alpha is the momentum compaction factor
   
   | See also :func:`atRFacc`

.. py:function:: atBunchLength(ring,ib,zn)

   | bunch length due to the potential well effect
   |  the output is the zerocurrent bunch length x bunch lengthening
   
   |    **bl = atBunchLength(ring,ib,zn)**
   
   |  Ib is the bunch current [A] (it may be a vector for multiple values)
   |  Zn is the longitudinal broadband impedance [Ohms]
   |  ring is the at ring without radiation
   |  BL is the bunch length in metres
   
   |    see also: BunchLength

.. py:function:: atRFacc(ring)

   | Computes RF acceptance of the ring
   |  **delta_max_rf = atRFacc(ring)**
   |    The functions computes the RF acceptance of the ring
   |    ring is tha at lattice without radiation
   |    delta_max_rf is the RF acceptance
   
   | See also :func:`RFacc`

.. py:function:: atSetCavityPhase

   | SETCAVITYPHASE     Set the TimeLag attribute of RF cavities
   
   | NEWRING=SETCAVITYPHASE(RING)
   |    Adjust the TimeLag attribute of RF cavities based on frequency,
   |    voltage and energy loss per turn, so that the synchronous phase is zero.
   |    An error occurs if all cavities do not have the same frequency.
   
   | NEWRING=SETCAVITYPHASE(...,'refpts',CAVPTS)
   |    CAVPTS is the location of RF cavities. This allows to ignore harmonic
   |    cavities.
   
   | WARNING: This function modifies the time reference,
   | this should be avoided
   
   | NEWRING=SETCAVITYPHASE(...,'method',METHOD)
   |    Choose the method for computing the energy loss per turn
   
   |  METHOD:   'integral': (default) The losses are obtained from
   |                        Losses = Cgamma / 2pi * EGeV^4 * I2
   |                        Takes into account bending magnets and wigglers.
   |            'tracking': The losses are obtained by tracking without cavities.
   |                        Needs radiation ON, takes into account all radiating elements.

.. py:function:: atsetcavity(ring,...,'frequency',frequency,...)

   | ATSECAVITY Set the parameters of RF cavities
   
   | **atsetcavity** may be used in two modes:
   
   | Upgrade mode
   | ===================================================
   |  By default, **atsetcavity** will act on the "main" cavities: they are defined by the
   |  cavpts ring property, or if absent by cavities at the lowest frequency.
   
   | **newring=atsetcavity(ring,...,'frequency',frequency,...)**
   |    Set the cavity frequency [Hz]. FREQUENCY is a scalar or an array as
   |    long as the list of selected cavities
   
   | **newring=atsetcavity(ring,...,'frequency','nominal',...)**
   |    Set the cavity frequency to the nominal value according to
   |    circumference and harmonic number
   
   | **newring=atsetcavity(ring,...,'frequency','nominal','dp',dp)**
   |    Set the cavity frequency to the nominal value for the specified dp
   
   | **newring=atsetcavity(ring,...,'frequency','nominal','dct',dct)**
   |    Set the cavity frequency to the nominal value for the specified dct
   
   | **newring=atsetcavity(ring,...,'frequency','nominal','df',df)**
   |    Set the cavity frequency to the nominal value + df
   
   | **newring=atsetcavity(ring,...,'voltage',voltage,...)**
   |    Set the total voltage (all cells) [V]. VOLTAGE will be distributed over the
   |    cells: CELL_VOLTAGE = VOLTAGE / PERIODICITY.
   |    Then if CELL_VOLTAGE is a scalar, it will be equally shared among the
   |    selected cavities. Otherwise it is an array as long as the list of
   |    selected cavities.
   
   | **newring=atsetcavity(ring,...,'harmnumber',h,...)**
   |    Set the harmonic number. H is a scalar or an array as
   |    long as the list of selected cavities
   
   | **newring=atsetcavity(ring,...,'timelag',timelag,...)**
   |    Set the time lag [m], . TIMELAG is a scalar or an array as
   |    long as the list of selected cavities
   
   | **newring=atsetcavity(ring,...,'cavpts',cavpts)**
   |    CAVPTS is the location of the selected RF cavities. The default is to act on the
   |    "main" cavities: they are defined by the cavpts ring property, or if absent by
   |    cavities at the lowest frequency.
   
   |   NOTES
   |   1. In this mode, the radiation state of the lattice is not modified.
   
   
   | Compatibility mode
   | ===================================================
   | **newring = atsetcavity(ring,rfv,radflag,harm_number)**
   |   RING         Ring structure
   |   RFV          RF voltage (full ring) [V]
   |   RADFLAG      0/1: activate/desactivate radiation (atradon/atradoff)
   |   HARMNUMBER 	Harmonic number (full ring)
   
   |   NOTES
   |   1. This mode is deprecated and should be replaced by
   |        **ring=atsetcavity(ring,'frequency','nominal','harmnumber',harm_number, 'voltage',rfv)**
   |        RING=atSetCavityPhase(RING) (optional)
   |        RING=atenable_6d(RING)      (optional)
   |   2. All the N cavities will have a voltage RFV/N
   |   3. sets the synchronous phase of the cavity assuming radiation is turned
   |      on radflag says whether or not we want radiation on, which affects
   |      synchronous phase.
   
   | See also :func:`atSetCavityPhase`, :func:`atsetRFCavity`, :func:`atenable_6d`, :func:`atdisable_6d`, :func:`atgetU0`

.. py:function:: cavityon(energy)

   | turns Cavities ON
   
   |  **cavityon** looks for elements that have field Frequency
   |     and sets PassMethod for them to RFCavityPass
   |  **cavityon(energy)**
   |     In addition sets the E0 field of the global variable GLOBVAL
   |     to energy - design energy [eV]
   |     If GLOBVAL does not exist **cavityon** creates it
   
   | See also :func:`cavityoff`, :func:`radiationon`, :func:`radiationoff`, :func:`setcavity`

.. py:function:: mcf(ring)

   | momentum compaction factor
   |  **mcf(ring)** calculates momentum compaction factor of RING
   
   |  **mcf(ring,dpp)** computes the momentum compaction for off-momentum DPP
   
   |  IMPORTANT!!!
   |  **mcf** gives a wrong result with 6-d rings. The RING should be set to 4d.
   | See also :func:`atdisable_6d`, :func:`check_6d`

.. py:function:: nus

   | Computes synchrotron tune from RF parameters
   |  **nus = nus** (VrfMV, alpha, U0MeV, E0MeV, h)
   |  this function return the synchrotron tune
   |  input:
   |  VrfMV is the RF voltage in MV
   |  alpha is the momentum compaction factor
   |  U0MeV is the energy lost per turn in MeV
   |  E0MeV is the beam energy in MeV
   |  h is the harmonic number

.. py:function:: phis

   | phase = phis(U0MeV,VrfMV)
   
   |  this function returns the synchronous phase in radians
   |  input:
   |  U0MeV is energy loss per turn in MeV
   |  VrfMV is the RF voltage in MV

