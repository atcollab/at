.. _mad82madx_module:

MAD82MADX
=========

.. rubric:: Functions


.. list-table::

   * - :func:`mad8TOmadx`
     - converts mad8 sequence files to madX

.. py:function:: mad8TOmadx(seqfilemad8)

   | converts mad8 sequence files to madX
   
   | function **[seqfilemadx]=mad8TOmadx(seqfilemad8)**
   
   |  This procedure reads a saved sequence in
   |  mad8 (SAVE,FILE='seqfilemad8';)
   |  and converts it to madx sequence
   |  every = goes to :=
   |  the order of the declarations is the same in the two files.
   
   |  works also with single mad8 files not containing comands, only
   |  definitions.
   |  does not translate call to files since those may change name
   
   |  parameters:
   |     - seqdilemad8=name of the mad8 lattice file
   |     - periodname (optional)= name of the period to use in madx (default is the filename)
   
   |  Simone Maria Liuzzo PhD@LNF 25-11-2011
   |      update 29-2-2012 : corrected a bug that would not translate correctly
   |      markers, kickers and monitor declared only by element name ("BPM: monitor" would not convet properly)

