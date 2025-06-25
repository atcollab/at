.. _madx2at_module:

MADX2AT
=======

.. rubric:: Functions


.. list-table::

   * - :func:`ParseAtributesMADX_2_AT`
     - determines atribute and sets field in sxs{i} structure AT
   * - :func:`atfrommadx`
     - function atfrommadx(seqfilemadX,E0,outfilename)
   * - :func:`buildATLattice`
     - given a list (cell array) of elements with specified field Spos (center of element (madx default)) in a
   * - :func:`reshapeToCellArray`
     - if CEL_CEL is a cell array of structures and cell arrays it converts it a

.. py:function:: ParseAtributesMADX_2_AT

   | determines atribute and sets field in sxs{i} structure AT
   
   |  created 6-sept-2012

.. py:function:: atfrommadx

   | function atfrommadx(seqfilemadX,E0,outfilename)
   |  tansform madX sequence file (savesequence) file into AT lattice structure.
   
   |  This procedure reads a saved lattice (sequence in madx) in madX
   |  and converts it to an AT lattice
   
   |  (madx comands to save the sequences :
   
   |   _______ MADX code _________
   |   use,period=sequencename1;
   |   use,period=sequencename2;
   |   use,period=sequencename2;
   |   SAVE,FILE='seqfilemadX.seq';
   |   ___________________________
   
   |   seqfilemadX.seq will contain sequencename1 sequencename2 sequencename3
   |   in the correct format in a single file
   
   |  )
   
   |  The routine outputs a Matlab macro with all the AT defitions and variables as
   |  in the madX file
   
   |  The order of the declarations is the same in the two files.
   |  declarations that contain other variables are moved to the end. (this may not be enough)
   
   
   |  Works also with single madX files not containing comands, only
   |  definitions.
   
   |  parameters:
   |     - seqfilemadX=name of the mad8 lattice file
   |     - E0  = design energy
   |     - outfilename (default: seqfilemadX_AT_LATTICE.mat)
   
   |  default pass methods:
   |           quadrupoles : StrMPoleSymplectic4Pass
   |           dipole : BndMPoleSymplectic4Pass
   |           multipole : StrMPoleSymplectic4Pass
   |           sextupole : StrMPoleSymplectic4Pass
   |           thinmultipole : ThinMPolePass
   |           correctors : ThinMPolePass
   |           cavity : DriftPass
   

.. py:function:: buildATLattice

   | given a list (cell array) of elements with specified field Spos (center of element (madx default)) in a
   |  vector returns a cell array with elements without Spos field and
   |  appropriate Drifts spaces between. Drifts of the same length have the same name.

.. py:function:: reshapeToCellArray

   | if CEL_CEL is a cell array of structures and cell arrays it converts it a
   |  cell array of structures.

