.. _at2madx_module:

AT2MADX
=======

.. rubric:: Functions


.. list-table::

   * - :func:`AT_2_madX`
     - function [elelat,defs,lines]=AT_2_madX(AT_ring,linename)

.. py:function:: AT_2_madX

   | function [elelat,defs,lines]=AT_2_madX(AT_ring,linename)
   |  this functions converts the AT lattice AT_ring in madX format.
   
   |  a MADX LINE is generated.
   
   |  file ['' linename '_lattice.madx'] is generated contiaining the lattice
   |  (elelat)  elements definitions (defs) and the LINE (lines).
   |  no other comands are introduced
   
   |  to test in MADX run (replace linename with apropriate name)
   |  madx < linename_lattice.madx
   
   |  to test with twiss and plot add following lines at the end of file:
   |  beam;
   |  use, period=linename;
   |  twiss;
   |  plot, haxis=s, vaxis1=betx,bety, vaxis2=dx;
   

