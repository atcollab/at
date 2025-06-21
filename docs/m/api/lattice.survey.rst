.. _survey_module:

survey
======

.. rubric:: Functions


.. list-table::

   * - :func:`atgeometry`
     - Computes the 2-D position of all elements (no vertical bend)
   * - :func:`atgeometry3`
     - Computes the 3-D position of all elements

.. py:function:: atgeometry(ring,refpts)

   | Computes the 2-D position of all elements (no vertical bend)
   | 
   | **posdata=atgeometry(ring,refpts)**
   | 
   | RING:	 AT structure
   | REFPTS: observation points (array of indexes or logical mask)
   |         The allowed range is 1 to length(RING)+1
   |         Defaults to 1:length(RING)+1
   | 
   | POSDATA:Structure array, same length as REFPTS, with 5 fields:
   |            x, y, angle, long, trans
   | 
   | **[posdata,radius]=atgeometry(ring,refpts)**
   |        Outputs the machine radius at the beginning of the lattice.
   |        Note: this is different from the radius usually defined as
   |        circumference/2/pi
   | 
   | **posdata=atgeometry(...,'centered')**
   |        The offset is set so that the origin is at the centre of the ring
   | 
   | **posdata=atgeometry(ring,refpts,offset)**
   |        Start at x=offset(1), y=offset(2). Ignored if 'centered' is set.
   |        A scalar offset value is equivalent to [0 OFFSET].
   | 
   | **posdata=atgeometry(...,'hangle',h_angle)**
   |        Set the initial trajectory angle
   | 
   | 
   | See also: ATGEOMETRY3

.. py:function:: atgeometry3(ring,refpts)

   | Computes the 3-D position of all elements
   | 
   | **posdata=atgeometry3(ring,refpts)**
   | 
   | RING:	 AT structure
   | REFPTS: observation points (array of indexes or logical mask)
   |         The allowed range is 1 to length(RING)+1
   |         Defaults to 1:length(RING)+1
   | 
   | POSDATA:Structure array, same length as REFPTS, with 3 fields:
   |            x, y, z
   | 
   | **posdata=atgeometry3(ring,refpts,offset)**
   |        Start at x=offset(1), y=offset(2), z=offset(3)
   |        a scalar offset value is equivalent to [0 OFFSET 0]
   | 
   | **posdata=atgeometry3(...,'hangle',h_angle)**
   |        Set the initial horizontal trajectory angle
   | 
   | **posdata=atgeometry3(...,'vangle',h_angle)**
   |        Set the initial vertical trajectory angle
   | 
   | See also: ATGEOMETRY

