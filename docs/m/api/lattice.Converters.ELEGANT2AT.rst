.. _elegant2at_module:

ELEGANT2AT
==========

.. py:module:: lattice.Converters.ELEGANT2AT

   ELEGANT2AT

.. rubric:: Functions


.. list-table::

   * - :func:`ParseAtributesELEGANT_2_AT`
     - determines atribute and sets field in sxs{i} structure AT
   * - :func:`elegant2at`
     - function (elegantlattice,E0,outfilename)

.. py:function:: ParseAtributesELEGANT_2_AT

   | determines atribute and sets field in sxs{i} structure AT
   
   |  created 6-sept-2012

.. py:function:: elegant2at

   | function (elegantlattice,E0,outfilename)
   |  tansform elegant %s.new file (save_lattice with output_seq=0) file into AT lattice structure.
   
   |  This procedure reads a saved elegant lattice and converts it to an AT lattice
   
   |  Elegant command to save the lattice sequence :
   
   |   _______ filename.ele_________
   |   &save_lattice
   |    filename = %s.new,
   |    output_seq=0,
   |  &end
   |   ___________________________
   
   |   filename.new will contain the list of all the elements that make up the
   |   lattice in the correct format in a single file
   
   |  The routine outputs a Matlab macro with all the AT defitions and variables as
   |  in the elegant file
   
   |  Works also with single elegant file not containing commands, only
   |  definitions.
   
   |  parameters:
   |     - elegantlattice = name of the elegant lattice file
   |     - E0  = design energy
   |     - outfilename (default: elegantlattice_AT_LATTICE.mat)
   
   |  default pass methods:
   |           quadrupoles : StrMPoleSymplectic4Pass
   |           dipole : BndMPoleSymplectic4Pass
   |           multipole : StrMPoleSymplectic4Pass
   |           sextupole : StrMPoleSymplectic4Pass
   |           thinmultipole : ThinMPolePass
   |           correctors : ThinMPolePass
   |           cavity : DriftPass
   

