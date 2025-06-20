.. _private_module:

private
=======

.. rubric:: Functions


.. list-table::

   * - :func:`decodeatargs`
     - DECODEATARGS separates arguments and resources

.. py:function:: decodeatargs

   | DECODEATARGS separates arguments and resources
   | 
   |   [RSRC,ARGS]=decodeatargs(DEFARGS,ARGLIST)
   | 
   |   INPUTS
   |     1. DEFARGS - Default values for mandatory argument
   |     2. ARGLIST - Arguments
   | 
   |   OUPUTS
   |     1. RSRC    - Optional arguments (all remaining arguments)
   |     2. ARGS    - Mandatory arguments
   | 
   |   See also getoption, getflag

