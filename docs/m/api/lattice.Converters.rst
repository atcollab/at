.. _converters_module:

Converters
==========

.. toctree::
   :hidden:

   lattice.Converters.AT2Elegant
   lattice.Converters.AT2G4BL
   lattice.Converters.AT2MAD8
   lattice.Converters.AT2MADX
   lattice.Converters.AT2OPA
   lattice.Converters.ELEGANT2AT
   lattice.Converters.MAD82MADX
   lattice.Converters.MADX2AT
   lattice.Converters.MADX2G4BL

.. rubric:: Modules


.. list-table::

   * - :ref:`at2elegant_module`
     - AT2ELEGANT
   * - :ref:`at2g4bl_module`
     - AT2G4BL
   * - :ref:`at2mad8_module`
     - AT2MAD8
   * - :ref:`at2madx_module`
     - AT2MADX
   * - :ref:`at2opa_module`
     - AT2OPA
   * - :ref:`elegant2at_module`
     - ELEGANT2AT
   * - :ref:`mad82madx_module`
     - MAD82MADX
   * - :ref:`madx2at_module`
     - MADX2AT
   * - :ref:`madx2g4bl_module`
     - MADX2G4BL

.. rubric:: Functions


.. list-table::

   * - :func:`readmad`
     - reads the file output of MAD commands

.. py:function:: readmad(filename)

   | reads the file output of MAD commands
   |  TWISS, STRUCTURE, SURVEY.
   | 
   |  **atlattice = readmad(filename)**
   | 
   |  **readmad** reads the MAD file header to determine the number of elements
   |  in the lattice, symmetry flag, the number of supperperiods etc.
   | 
   |  Then it interprets the entry for each element in the MAD output file.
   |  The topology of the lattice is completely determined by
   |  Length, Bending Angle, and Ttilt Angle in each element
   | 
   |  **readmad** uses MAD TYPES and the values of to determine
   |  which pass-method function in AT to use.
   | 
   |  MAD TYPE      |  AT PassMethod
   |  ----------------------------------
   |  DRIFT         |  DriftPass
   |  SBEND         |  BendLinearPass, BndMPoleSymplectic4Pass
   |  QUADRUPOLE    |  QualdLinearPass
   |  SEXTUPOLE     |  StrMPoleSymplectic4Pass
   |  OCTUPOLE      |  StrMPoleSymplectic4Pass
   |  MULTIPOLE     |  !!! Not implemented, in future - ThinMPolePass
   |  RFCAVITY      |  RFCavityPass
   |  KICKER        |  CorrectorPass
   |  HKICKER       |  CorrectorPass
   |  VKICKER       |  CorrectorPass
   |  MONITOR       |  IdentityPass
   |  HMONITOR      |  IdentityPass
   |  VMONITOR      |  IdentityPass
   |  MARKER        |  IdentityPass
   |  -----------------------------------
   |  all others    |  Length=0 -> IdentityPass, Length~=0 -> DriftPass

