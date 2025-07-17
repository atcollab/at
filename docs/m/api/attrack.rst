.. _attrack_module:

attrack
=======

.. py:module:: attrack

   Tracking functions

.. rubric:: Functions


.. list-table::

   * - :func:`atpass`
     - Low-level tracking engine
   * - :func:`elempass`
     - Tracks particles through a single element
   * - :func:`gpuinfo`
     - Get information on GPUs available for tracking
   * - :func:`gpupass`
     - Low-level tracking engine on GPU
   * - :func:`linepass`
     - tracks particles through each element of the cell array LINE
   * - :func:`ringpass`
     - tracks particles through each element of the cell array RING

.. py:function:: atpass(lattice,rin,mode,nturns,refpts,prefunc,postfunc,nhist,numthreads,ringprops)

   | Low-level tracking engine
   
   |  **rout = atpass(lattice,rin,mode,nturns,refpts,prefunc,postfunc,nhist,numthreads,ringprops)**
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
   
   |  **[rout,lost] = atpass(...)**
   |    Additionally return information on lost particles
   
   |    LOST        1x1 structure with the following fields:
   |                turn        1xN vector, turn number where the particle is lost
   |                element     1xN vector, element number where the particle is lost
   |                coordinates 6xNxNHIST matrix, coordinates at the entrance of the
   |                            element where the particle was lost
   
   | See also :func:`ringpass`, :func:`linepass`

.. py:function:: elempass(elem,rin)

   | Tracks particles through a single element
   
   | **rout=elempass(elem,rin)** Tracks particles through ELEM
   
   |    ELEM:       lattice element
   |    RIN         6xN matrix: input coordinates of N particles
   
   |    ROUT        6xN matrix: output coordinates of N particles
   
   |   **rout=elempass(...,'passmethod',passmethod,...)**
   |      Use PASSMETHOD (default: ELEM.PassMethod)
   
   |   **rout=elempass(...,'energy',energy,...)**
   |      Use ENERGY and ignore the 'Energy' field of elements
   
   |   **rout=elempass(...,'particle',particle,...)**
   |      Use PARTICLE (default: relativistic)
   
   | See also :func:`ringpass`, :func:`linepass`

.. py:function:: gpuinfo

   | Get information on GPUs available for tracking
   
   |  **info = gpuinfo**
   
   |  INFO:    1xn structure with the following fields:
   |           Name:       GPU name
   |           Version:    CUDA compute capability (? for OpenCL)
   |           CoreNumber: Multi processor number
   |           Platform:   Platform name

.. py:function:: gpupass(lattice,rin,mode,nturns,refpts,turn,keepcounter,gpupool,integrator)

   | Low-level tracking engine on GPU
   
   |  **rout = gpupass(lattice,rin,mode,nturns,refpts,turn,keepcounter,gpupool,integrator)**
   |    LATTICE     AT lattice
   |    RIN         6xN matrix: input coordinates of N particles
   |    MODE        0 - reuse lattice
   |                1 - new lattice
   |    NTURNS      number of turns
   |    REFPTS      Indexes of elements where the trajectory is observed
   |                May run from 1 to length(LATTICE)+1
   |    KEEPCOUNTER 0 - Start counting turn from 0
   |                1 - Keep last turn counter of the previous **gpupass** call
   |    GPUPOOL     GPU to use (see gpuinfo)
   |    INTEGRATOR  Type of integrator to use
   |                1: Euler 1st order, 1 drift/1 kick per step
   |                2: Verlet 2nd order, 1 drift/2 kicks per step
   |                3: Ruth 3rd order, 3 drifts/3 kicks per step
   |                4: Forest/Ruth 4th order, 4 drifts/3 kicks per step (Default)
   |                5: Optimal 4th order from R. Mclachlan, 4 drifts/4 kicks per step
   |                6: Yoshida 6th order, 8 drifts/7 kicks per step
   
   |  **[rout,lost] = gpupass(...)**
   |    Additionally return information on lost particles
   
   |    LOST        1x1 structure with the following fields:
   |                turn                1xN vector, turn number where the particle is lost
   |                element             1xN vector, element number where the particle is lost
   |                coordinates_at_loss 6xN matrix, coordinates at the entrance of the
   |                                    element where the particle was lost
   
   | See also :func:`ringpass`, :func:`linepass`

.. py:function:: linepass(line,rin) tracks particle(s)

   | tracks particles through each element of the cell array LINE
   |  calling the element-specific tracking function specified in the
   |  LINE{i}.PassMethod field.
   
   |  **rout=linepass(line,rin) tracks particle(s)** with initial
   |     condition(s) RIN for NTURNS turns to the end of the LINE
   
   |    LINE        AT lattice
   |    RIN         6xN matrix: input coordinates of N particles
   
   |    ROUT        6xN matrix: output coordinates of N particles at
   |                the end of LINE
   
   |  **rout=linepass(line,rin,refpts)** also returns intermediate results
   |      at the entrance of each element specified in the REFPTS
   
   |     REFPTS is an array of increasing indexes that selects elements
   |      between 1 and length(LINE)+1.
   |      See further explanation of REFPTS in the 'help' for FINDSPOS
   |    ROUT        6x(N*length(REFPTS)) matrix: output coordinates of N particles at
   |                each reference point
   
   |      NOTE:
   |      **linepass(line,rin,length(line)+1) is the same as linepass(line,rin)**
   |      since the reference point length(LINE)+1 is the exit of the last element
   |      **linepass(line,rin,1)** is a copy of RIN since the
   |      reference point 1 is the entrance of the first element
   
   |  **[rout, lost]=linepass(...)**
   |   Return additionally an information on lost particles
   |     LOST	1xN logical vector, indicating lost particles
   |     If only one output is given, loss information is saved in
   |     global variable LOSSFLAG
   
   |  **[rout, loss, lossinfo]=linepass(...,'nhist',nhist,...)**
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
   
   |  **rout=linepass(...,'keeplattice')** Tracking with the 'KeepLattice' flag is
   |    more efficient because it reuses persistent data structures stored in
   |    memory in previous calls to **linepass**.
   
   | 	!!! In order to use this option, **linepass** must first be called
   | 	without the 'KeepLattice' flag. It then assumes that the elements in LINE
   |  	DO NOT CHANGE between calls. Otherwise, **linepass** must be called again
   |    without 'KeepLattice'.
   
   |  **rout=linepass(...,'reuse')** is kept for compatibilty with previous
   |  versions. It has no effect.
   
   |  **rout=linepass(...,'seed',seed)**  The random generators are reset to start
   |    with SEED.
   
   |  **rout=linepass(...,'omp_num_threads',nthreads)**  Number of OpenMP threads.
   |    By default, OpenMP chooses the number of threads.
   
   |  **rfin=linepass(...,prefunc)**
   |  **rfin=linepass(...,prefunc,postfunc)**
   |  **rfin=linepass(...,cell(0),postfunc)**
   |     PREFUNC and POSTFUNC are function handles, PREFUNC is called
   |     immediately before tracking each element, POSTFUNC is called
   |     immediately after each element. Functions are called as:
   
   |        ROUT=FUNC(ELEMENT, RIN, NTURN, NELEMENT)
   
   |    and is allowed to modify the particle coordinates
   
   | See also :func:`ringpass`

.. py:function:: ringpass(ring,rin,nturns) tracks particle(s)

   | tracks particles through each element of the cell array RING
   |  calling the element-specific tracking function specified in the
   |  RING{i}.PassMethod field.
   
   |  **rout=ringpass(ring,rin,nturns) tracks particle(s)** with initial
   |     condition(s) RIN for NTURNS turns
   
   |    RING        AT lattice
   |    RIN         6xN matrix: input coordinates of N particles
   |    NTURNS      Number of turns to perform (default: 1)
   
   |    ROUT        6x(N*NTURNS) matrix: output coordinates of N particles at
   |                the exit of each turn
   
   |  **[rout, lost]=ringpass(...)**
   |   Return additionally an information on lost particles
   |     LOST	1xN logical vector, indicating lost particles
   |     If only one output is given, loss information is saved in
   |     global variable LOSSFLAG
   
   |  **[rout, lost, nturns]=ringpass(...)**
   |   Return additionally the number of turns performed by each particle
   | 	NTURNS	1xN vector, number of turns performed
   
   |  **[rout, loss, nturns, lossinfo]=ringpass(...,'nhist',nhist,...)**
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
   
   |  **rout=ringpass(...,'keeplattice')** Tracking with the 'KeepLattice' flag is
   |    more efficient because it reuses persistent data structures stored in
   |    memory in previous calls to **ringpass**.
   
   | 	!!! In order to use this option, **ringpass** must first be called
   | 	without the 'KeepLattice' flag. It then assumes that the elements in RING
   | 	DO NOT CHANGE between calls. Otherwise, **ringpass**  must be called again
   |    without 'KeepLattice'.
   
   |  **rout=ringpass(...,'reuse')** is kept for compatibilty with previous
   |  versions. It has no effect.
   
   |  **rout=ringpass(...,'seed',seed)**  The random generators are reset to start
   |    with SEED.
   
   |  **rout=ringpass(...,'turn',turn)**    Initial turn number. Default 0.
   |    The turn number is necessary to compute the absolute path length used
   |    by RFCavityPass. Ignored if KeepCounter is set.
   
   |  **rout=ringpass(...,'keepcounter')**  The turn number starts with the last
   |    turn of the previous call.
   
   |  NOTE:
   |  To resume an interrupted tracking (for instance to get intermediate
   |  results), one must use one of the 'turn' option or 'KeepCounter' flag to
   |  ensure the continuity of the turn number.
   
   |  **rout=ringpass(...,'omp_num_threads',nthreads)**  Number of OpenMP threads.
   |    By default, OpenMP chooses the number of threads.
   
   |  **rout=ringpass(...,'silent')** does not output the particle coordinates at
   |     each turn but only at the end of the tracking
   
   |  **rout=ringpass(...,prefunc)**
   |  **rout=ringpass(...,prefunc,postfunc)**
   |  **rout=ringpass(...,cell(0),postfunc)**
   |    PREFUNC and POSTFUNC are function handles, PREFUNC is called
   |    immediately before tracking each element, POSTFUNC is called
   |    immediately after each element. Functions are called as:
   
   |        ROUT=FUNC(ELEMENT, RIN, NTURN, NELEMENT)
   
   |    and are allowed to modify the particle coordinates
   
   | See also :func:`linepass`

