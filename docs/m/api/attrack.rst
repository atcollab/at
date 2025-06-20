.. _attrack_module:

attrack
=======

.. rubric:: Functions


.. list-table::

   * - :func:`gpuinfo`
     -  INFO = at_gpuinfo()
   * - :func:`gpupass`
     - ATPASS is a numerical tracking engine for AT
   * - :func:`atpass`
     - ATPASS is a numerical tracking engine for AT
   * - :func:`elempass`
     - ELEMPASS Tracks particles through a single element
   * - :func:`ringpass`
     - RINGPASS tracks particles through each element of the cell array RING
   * - :func:`linepass`
     - LINEPASS tracks particles through each element of the cell array LINE

.. py:function:: gpuinfo

   |  INFO = at_gpuinfo()
   |  INFO          1xn structure with the following fields:
   |                Name: GPU name
   |                Version: CUDA compute capability (? for OpenCL)
   |                CoreNumber: Multi processor number
   |                Platform: Platform name

.. py:function:: gpupass

   | ATPASS is a numerical tracking engine for AT
   | 
   |  ROUT = GPUPASS(LATTICE,RIN,MODE,NTURNS,REFPTS,TURN,KEEPCOUNTER,GPUPOOL,INTEGRATOR)
   |    LATTICE     AT lattice
   |    RIN         6xN matrix: input coordinates of N particles
   |    MODE        0 - reuse lattice
   |                1 - new lattice
   |    NTURNS      number of turns
   |    REFPTS      Indexes of elements where the trajectory is observed
   |                May run from 1 to length(LATTICE)+1
   |    KEEPCOUNTER 0 - Start counting turn from 0
   |                1 - Keep last turn counter of the previous gpupass call
   |    GPUPOOL     GPU to use (see gpuinfo)
   |    INTEGRATOR  Type of integrator to use
   |                1: Euler 1st order, 1 drift/1 kick per step
   |                2: Verlet 2nd order, 1 drift/2 kicks per step
   |                3: Ruth 3rd order, 3 drifts/3 kicks per step
   |                4: Forest/Ruth 4th order, 4 drifts/3 kicks per step (Default)
   |                5: Optimal 4th order from R. Mclachlan, 4 drifts/4 kicks per step
   |                6: Yoshida 6th order, 8 drifts/7 kicks per step
   | 
   |  [ROUT,LOST] = GPUPASS(...)
   |    Additionally return information on lost particles
   | 
   |    LOST        1x1 structure with the following fields:
   |                turn                1xN vector, turn number where the particle is lost
   |                element             1xN vector, element number where the particle is lost
   |                coordinates_at_loss 6xN matrix, coordinates at the entrance of the
   |                                    element where the particle was lost
   | 
   |  See also  RINGPASS, LINEPASS

.. py:function:: atpass

   | ATPASS is a numerical tracking engine for AT
   | 
   |  ROUT = ATPASS(LATTICE,RIN,MODE,NTURNS,REFPTS,PREFUNC,POSTFUNC,NHIST,NUMTHREADS,RINGPROPS)
   |    LATTICE     AT lattice
   |    RIN         6xN matrix: input coordinates of N particles
   |    MODE        0 - reuse lattice
   |                1 - new lattice
   |    NTURNS      number of turns
   |    REFPTS      Indexes of elements where the trajectory is observed
   |                May run from 1 to length(LATTICE)+1
   |    PREFUNC     function called before each element
   |    POSTFUNC    function called after each element
   |    NHIST       length of history buffer. Optional, default 1
   |    NUMTHREADS  Number of threads in OpenMP. Optional, default: automatic
   |    RINGPROPS   Ring properties (energy and particle).
   | 
   |  [ROUT,LOST] = ATPASS(...)
   |    Additionally return information on lost particles
   | 
   |    LOST        1x1 structure with the following fields:
   |                turn        1xN vector, turn number where the particle is lost
   |                element     1xN vector, element number where the particle is lost
   |                coordinates 6xNxNHIST matrix, coordinates at the entrance of the
   |                            element where the particle was lost
   | 
   |  See also  RINGPASS, LINEPASS

.. py:function:: elempass

   | ELEMPASS Tracks particles through a single element
   | 
   | ROUT=ELEMPASS(ELEM,RIN) Tracks particles through ELEM
   | 
   |    ELEM:       lattice element
   |    RIN         6xN matrix: input coordinates of N particles
   | 
   |    ROUT        6xN matrix: output coordinates of N particles
   | 
   |   ROUT=ELEMPASS(...,'PassMethod',PASSMETHOD,...)
   |      Use PASSMETHOD (default: ELEM.PassMethod)
   | 
   |   ROUT=ELEMPASS(...,'Energy',ENERGY,...)
   |      Use ENERGY and ignore the 'Energy' field of elements
   | 
   |   ROUT=ELEMPASS(...,'Particle',PARTICLE,...)
   |      Use PARTICLE (default: relativistic)
   | 
   |  See also: RINGPASS, LINEPASS

.. py:function:: ringpass

   | RINGPASS tracks particles through each element of the cell array RING
   |  calling the element-specific tracking function specified in the
   |  RING{i}.PassMethod field.
   | 
   |  ROUT=RINGPASS(RING,RIN,NTURNS) tracks particle(s) with initial
   |     condition(s) RIN for NTURNS turns
   | 
   |    RING        AT lattice
   |    RIN         6xN matrix: input coordinates of N particles
   |    NTURNS      Number of turns to perform (default: 1)
   | 
   |    ROUT        6x(N*NTURNS) matrix: output coordinates of N particles at
   |                the exit of each turn
   | 
   |  [ROUT, LOST]=RINGPASS(...)
   |   Return additionally an information on lost particles
   |     LOST	1xN logical vector, indicating lost particles
   |     If only one output is given, loss information is saved in
   |     global variable LOSSFLAG
   | 
   |  [ROUT, LOST, NTURNS]=RINGPASS(...)
   |   Return additionally the number of turns performed by each particle
   | 	NTURNS	1xN vector, number of turns performed
   | 
   |  [ROUT, LOSS, NTURNS, LOSSINFO]=RINGPASS(...,'nhist',NHIST,...)
   |   Return additional information on lost particles
   |    NHIST       number elements before the loss to be traced (default: 1)
   |    LOSSINFO	1x1 structure with the following fields:
   |                lost                 1xN logical vector, indicating lost particles
   |                turn                 1xN vector, turn number where the particle is lost
   |                element              1xN vector, element number where the particle is lost
   |                coordinates_at_loss  6xN array, coordinates at the exit of
   |                                     the element where the particle is lost
   |                                     (sixth coordinate is inf if particle is lost in a physical aperture)
   |                coordinates          6xNxNHIST array, coordinates at the entrance of the
   |                                     LHIST elements before the loss
   | 
   |  ROUT=RINGPASS(...,'KeepLattice') Tracking with the 'KeepLattice' flag is
   |    more efficient because it reuses persistent data structures stored in
   |    memory in previous calls to RINGPASS.
   | 
   | 	!!! In order to use this option, RINGPASS must first be called
   | 	without the 'KeepLattice' flag. It then assumes that the elements in RING
   | 	DO NOT CHANGE between calls. Otherwise, RINGPASS  must be called again
   |    without 'KeepLattice'.
   | 
   |  ROUT=RINGPASS(...,'reuse') is kept for compatibilty with previous
   |  versions. It has no effect.
   | 
   |  ROUT=RINGPASS(...,'seed',SEED)  The random generators are reset to start
   |    with SEED.
   | 
   |  ROUT=RINGPASS(...,'turn',TURN)    Initial turn number. Default 0.
   |    The turn number is necessary to compute the absolute path length used
   |    by RFCavityPass. Ignored if KeepCounter is set.
   | 
   |  ROUT=RINGPASS(...,'KeepCounter')  The turn number starts with the last
   |    turn of the previous call.
   | 
   |  NOTE:
   |  To resume an interrupted tracking (for instance to get intermediate
   |  results), one must use one of the 'turn' option or 'KeepCounter' flag to
   |  ensure the continuity of the turn number.
   | 
   |  ROUT=RINGPASS(...,'omp_num_threads',NTHREADS)  Number of OpenMP threads.
   |    By default, OpenMP chooses the number of threads.
   | 
   |  ROUT=RINGPASS(...,'Silent') does not output the particle coordinates at
   |     each turn but only at the end of the tracking
   | 
   |  ROUT=RINGPASS(...,PREFUNC)
   |  ROUT=RINGPASS(...,PREFUNC,POSTFUNC)
   |  ROUT=RINGPASS(...,cell(0),POSTFUNC)
   |    PREFUNC and POSTFUNC are function handles, PREFUNC is called
   |    immediately before tracking each element, POSTFUNC is called
   |    immediately after each element. Functions are called as:
   | 
   |        ROUT=FUNC(ELEMENT, RIN, NTURN, NELEMENT)
   | 
   |    and are allowed to modify the particle coordinates
   | 
   |  See also: LINEPASS

.. py:function:: linepass

   | LINEPASS tracks particles through each element of the cell array LINE
   |  calling the element-specific tracking function specified in the
   |  LINE{i}.PassMethod field.
   | 
   |  ROUT=LINEPASS(LINE,RIN) tracks particle(s) with initial
   |     condition(s) RIN for NTURNS turns to the end of the LINE
   | 
   |    LINE        AT lattice
   |    RIN         6xN matrix: input coordinates of N particles
   | 
   |    ROUT        6xN matrix: output coordinates of N particles at
   |                the end of LINE
   | 
   |  ROUT=LINEPASS(LINE,RIN,REFPTS) also returns intermediate results
   |      at the entrance of each element specified in the REFPTS
   | 
   |     REFPTS is an array of increasing indexes that selects elements
   |      between 1 and length(LINE)+1.
   |      See further explanation of REFPTS in the 'help' for FINDSPOS
   |    ROUT        6x(N*length(REFPTS)) matrix: output coordinates of N particles at
   |                each reference point
   | 
   |      NOTE:
   |      LINEPASS(LINE,RIN,length(LINE)+1) is the same as LINEPASS(LINE,RIN)
   |      since the reference point length(LINE)+1 is the exit of the last element
   |      LINEPASS(LINE,RIN,1) is a copy of RIN since the
   |      reference point 1 is the entrance of the first element
   | 
   |  [ROUT, LOST]=LINEPASS(...)
   |   Return additionally an information on lost particles
   |     LOST	1xN logical vector, indicating lost particles
   |     If only one output is given, loss information is saved in
   |     global variable LOSSFLAG
   | 
   |  [ROUT, LOSS, LOSSINFO]=LINEPASS(...,'nhist',NHIST,...)
   |   Return additional information on lost particles
   |    NHIST       number elements before the loss to be traced (default: 1)
   |    LOSSINFO	1x1 structure with the following fields:
   |                lost                 1xN logical vector, indicating lost particles
   |                turn                 1xN vector, turn number where the particle is lost
   |                element              1xN vector, element number where the particle is lost
   |                coordinates_at_loss  6xN array, coordinates at the exit of
   |                                     the element where the particle is lost
   |                                     (sixth coordinate is inf if particle is lost in a physical aperture)
   |                coordinates          6xNxNHIST array, coordinates at the entrance of the
   |                                     LHIST elements before the loss
   | 
   |  ROUT=LINEPASS(...,'KeepLattice') Tracking with the 'KeepLattice' flag is
   |    more efficient because it reuses persistent data structures stored in
   |    memory in previous calls to LINEPASS.
   | 
   | 	!!! In order to use this option, LINEPASS must first be called
   | 	without the 'KeepLattice' flag. It then assumes that the elements in LINE
   |  	DO NOT CHANGE between calls. Otherwise, LINEPASS must be called again
   |    without 'KeepLattice'.
   | 
   |  ROUT=LINEPASS(...,'reuse') is kept for compatibilty with previous
   |  versions. It has no effect.
   | 
   |  ROUT=LINEPASS(...,'seed',SEED)  The random generators are reset to start
   |    with SEED.
   | 
   |  ROUT=LINEPASS(...,'omp_num_threads',NTHREADS)  Number of OpenMP threads.
   |    By default, OpenMP chooses the number of threads.
   | 
   |  Rfin=LINEPASS(...,PREFUNC)
   |  Rfin=LINEPASS(...,PREFUNC,POSTFUNC)
   |  Rfin=LINEPASS(...,cell(0),POSTFUNC)
   |     PREFUNC and POSTFUNC are function handles, PREFUNC is called
   |     immediately before tracking each element, POSTFUNC is called
   |     immediately after each element. Functions are called as:
   | 
   |        ROUT=FUNC(ELEMENT, RIN, NTURN, NELEMENT)
   | 
   |    and is allowed to modify the particle coordinates
   | 
   |  See also: RINGPASS

