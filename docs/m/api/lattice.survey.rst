.. _survey_module:

survey
======

.. rubric:: Functions


.. list-table::

   * - :func:`atgeometry3`
     - ATGEOMETRY3 Computes the 3-D position of all elements
   * - :func:`atgeometry`
     - ATGEOMETRY Computes the 2-D position of all elements (no vertical bend)

.. py:function:: atgeometry3

   | ATGEOMETRY3 Computes the 3-D position of all elements
   | 
   | POSDATA=ATGEOMETRY3(RING,REFPTS)
   | 
   | RING:	 AT structure
   | REFPTS: observation points (array of indexes or logical mask)
   |         The allowed range is 1 to length(RING)+1
   |         Defaults to 1:length(RING)+1
   | 
   | POSDATA:Structure array, same length as REFPTS, with 3 fields:
   |            x, y, z
   | 
   | POSDATA=ATGEOMETRY3(RING,REFPTS,OFFSET)
   |        Start at x=offset(1), y=offset(2), z=offset(3)
   |        a scalar offset value is equivalent to [0 OFFSET 0]
   | 
   | POSDATA=ATGEOMETRY3(...,'Hangle',h_angle)
   |        Set the initial horizontal trajectory angle
   | 
   | POSDATA=ATGEOMETRY3(...,'Vangle',h_angle)
   |        Set the initial vertical trajectory angle
   | 
   | See also: ATGEOMETRY

.. py:function:: atgeometry

   | ATGEOMETRY Computes the 2-D position of all elements (no vertical bend)
   | 
   | POSDATA=ATGEOMETRY(RING,REFPTS)
   | 
   | RING:	 AT structure
   | REFPTS: observation points (array of indexes or logical mask)
   |         The allowed range is 1 to length(RING)+1
   |         Defaults to 1:length(RING)+1
   | 
   | POSDATA:Structure array, same length as REFPTS, with 5 fields:
   |            x, y, angle, long, trans
   | 
   | [POSDATA,RADIUS]=ATGEOMETRY(RING,REFPTS)
   |        Outputs the machine radius at the beginning of the lattice.
   |        Note: this is different from the radius usually defined as
   |        circumference/2/pi
   | 
   | POSDATA=ATGEOMETRY(...,'centered')
   |        The offset is set so that the origin is at the centre of the ring
   | 
   | POSDATA=ATGEOMETRY(RING,REFPTS,OFFSET)
   |        Start at x=offset(1), y=offset(2). Ignored if 'centered' is set.
   |        A scalar offset value is equivalent to [0 OFFSET].
   | 
   | POSDATA=ATGEOMETRY(...,'Hangle',h_angle)
   |        Set the initial trajectory angle
   | 
   | 
   | See also: ATGEOMETRY3

