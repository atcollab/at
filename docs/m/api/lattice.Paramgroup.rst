.. _paramgroup_module:

Paramgroup
==========

.. rubric:: Functions


.. list-table::

   * - :func:`mkparamgroup`
     - simplifies creation of AT parameter groups
   * - :func:`restoreparamgroup`
     - restores the values of multiple physical
   * - :func:`saveparamgroup`
     - saves the values of multiple physical
   * - :func:`setparamgroup`
     - modifies a group of parameters

.. py:function:: mkparamgroup(lattice,elemindex,paramstr)

   | simplifies creation of AT parameter groups
   |  It group one or more elements in the
   |  same family and simultaneously vary
   | 
   |  **mkparamgroup(lattice,elemindex,paramstr)**
   |  **mkparamgroup(lattice,famnamestr,paramstr)**
   |  **mkparamgroup(lattice,famnamestr,kidnum,paramstr)**
   | 
   |  LATTICE
   |  FAMNAMESTR
   | 
   | 
   |  PARAMSTR: 'TILT','K1','K2','K3'
   |  wjc 2-09-04 changed index 'i' to 'k'

.. py:function:: restoreparamgroup(lattice,paramgroup)

   | restores the values of multiple physical
   |  parameters of the lattice.
   |  **newlattice = restoreparamgroup(lattice,paramgroup)**
   | 
   |  See also: ATPARAMGROUP RESTORPARAMGROUP SAVEPARAMGROUP

.. py:function:: saveparamgroup(lattice,paramgroup)

   | saves the values of multiple physical
   |  parameters of the lattice in the special SavedValue field of
   |  AT parameter group structure. The result can be late used
   |  with RESTOREPARAMGROUP
   | 
   |  **paramgroup = saveparamgroup(lattice,paramgroup)**
   | 
   |  See also: ATPARAMGROUP RESTORPARAMGROUP SETPARAMGROUP

.. py:function:: setparamgroup(lattice,paramgroup,pvalue)

   | modifies a group of parameters
   |  **newlattice = setparamgroup(lattice,paramgroup,pvalue)**
   | 
   |  See also: ATPARAMGROUP RESTORPARAMGROUP

