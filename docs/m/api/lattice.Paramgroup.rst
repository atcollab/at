.. _paramgroup_module:

Paramgroup
==========

.. rubric:: Functions


.. list-table::

   * - :func:`saveparamgroup`
     - SAVEPARAMGROUP saves the values of multiple physical
   * - :func:`restoreparamgroup`
     - RESTOREPARAMGROUP restores the values of multiple physical
   * - :func:`setparamgroup`
     - SETPARAMGROUP modifies a group of parameters
   * - :func:`mkparamgroup`
     - MKPARAMGROUP simplifies creation of AT parameter groups

.. py:function:: saveparamgroup

   | SAVEPARAMGROUP saves the values of multiple physical
   |  parameters of the lattice in the special SavedValue field of
   |  AT parameter group structure. The result can be late used
   |  with RESTOREPARAMGROUP
   | 
   |  PARAMGROUP = saveparamgroup(LATTICE,PARAMGROUP)
   | 
   |  See also: ATPARAMGROUP RESTORPARAMGROUP SETPARAMGROUP

.. py:function:: restoreparamgroup

   | RESTOREPARAMGROUP restores the values of multiple physical
   |  parameters of the lattice.
   |  NEWLATTICE = RESTOREPARAMGROUP(LATTICE,PARAMGROUP)
   | 
   |  See also: ATPARAMGROUP RESTORPARAMGROUP SAVEPARAMGROUP

.. py:function:: setparamgroup

   | SETPARAMGROUP modifies a group of parameters
   |  NEWLATTICE = setparamgroup(LATTICE,PARAMGROUP,PVALUE)
   | 
   |  See also: ATPARAMGROUP RESTORPARAMGROUP

.. py:function:: mkparamgroup

   | MKPARAMGROUP simplifies creation of AT parameter groups
   |  It group one or more elements in the
   |  same family and simultaneously vary
   | 
   |  MKPARAMGROUP(LATTICE,ELEMINDEX,PARAMSTR)
   |  MKPARAMGROUP(LATTICE,FAMNAMESTR,PARAMSTR)
   |  MKPARAMGROUP(LATTICE,FAMNAMESTR,KIDNUM,PARAMSTR)
   | 
   |  LATTICE
   |  FAMNAMESTR
   | 
   | 
   |  PARAMSTR: 'TILT','K1','K2','K3'
   |  wjc 2-09-04 changed index 'i' to 'k'

