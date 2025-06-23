.. _touschekpiwinski_module:

TouschekPiwinski
================

.. rubric:: Functions


.. list-table::

   * - :func:`MomAperture_Project2Start`
     - calculates the local momentum aperture.
   * - :func:`MomAperture_allRing`
     - MomAperture_allRing returns positive and negative momentum aperture
   * - :func:`TLT_IntPiw`
     - integral in Piwinski Formula for the Lifetime
   * - :func:`TLT_IntPiw_k`
     - integral in Piwinski Formula for the Lifetime with u=tan^2(k)
   * - :func:`TouschekPiwinskiLifeTime`
     - function [Tl,contributionsTL]=TouschekPiwinskiLifeTime(ring,dpp,Ib,...)
   * - :func:`momentum_aperture_at`
     - 

.. py:function:: MomAperture_Project2Start(thering)

   | calculates the local momentum aperture.
   | 
   |  **momaperture_project2start** is a Bipartition search of the negative and
   |  positive stability thesholds in the 5th dimension (relative momentum).
   |   -The 6D closed orbit is taken into account.
   |   -Particles launched at different REFPTS along the ring are first projected
   |   to the ring last element so that all particles can be tracked together.
   | 
   |  **[etn, etp] = momaperture_project2start(thering)**
   | 
   |  Inputs:
   |        THERING: ring used for tracking.
   |  Options:
   |        REFPTS: REFPTS where to calculate the momentum acceptance.
   |                Default 1:numel(THERING);
   |        nturns: Number of turns to track. Default 1000
   |        dptol:  resolution in momentum acceptance. Default 1e-4
   |        dpuguess: unstable momentum threshold guess. Default [].
   |                If not given it uses the linear momentum acceptance delta_max
   |                from ringpara.
   |        troffset: [x y] starting transverse offset for the tracking.
   |                Default [1e-6 1e-6]
   |        verbose: boolean indicating verbose mode. Default false.
   |        epsilon6D: if not passed, all particles are tracked.
   |                If epsilon6D is given, we track for nturns only
   |                particles having 6D coordinates different by epsilon6D
   |                after being projected to the end of the ring.
   |  Output:
   |        ETN: stability threshold for positive off momentum particles
   |        ETP: stability threshold for negative off momentum particles
   | 
   | 
   |  Other functions in the file:
   | 
   |  Loste = Multiorigin_ringpass_islost
   |  Returns a boolean array: tells whether the particle launched at
   |  the reference point refpts with positive and negative momentum offset is
   |  lost or not.

.. py:function:: MomAperture_allRing(..., nturns)

   | MomAperture_allRing returns positive and negative momentum aperture
   |                     boundaries where the particle is still alive.
   | 
   |                     The boundary width (i.e. the uncertainty) is equal to
   |                     energystep / (splitdivisor^ntimessplit),
   |                     meaning that one more step of this size makes the
   |                     particle unstable.
   | 
   |  [map_l,map_h] ...
   |     = **momaperture_allring**(
   |                           THERING, ...
   |                           POINTS ...
   |                          )
   |  **[...] = momaperture_allring(..., nturns)**
   |          Tracks over NTURNS to get the momentum aperture. Default 100
   |          e.g **[dppm,dppp]=momaperture_allring(thering,positions,nturns)**
   | 
   |  **[...] = momaperture_allring(..., 'reforbit',orbitin)**
   |          The initial particle coordinates are taken from ORBITIN.
   |          Default zeros(6,length(POINTS))
   | 
   |  **[...] = momaperture_allring(..., 'xyinitoffsets',[x y])**
   |          The transverse offsets to add to the reference orbit.
   |          Default 1e-5*ones(length(POINTS),2)
   | 
   |  **[...] = momaperture_allring(..., 'deltalimits',[deltapos deltaneg])**
   |          The energy offset limits. Default [0.1 -0.1]
   | 
   |  **[...] = momaperture_allring(..., 'initialguess',[posguess negguess])**
   |          The starting point of the recursive energy offsets.
   |          Default [0.0 0.0]
   | 
   |  **[...] = momaperture_allring(..., 'energysteps',[posstep negstep])**
   |          The positive and negative initial energy steps.
   |          Default [0.01 -0.01]
   | 
   |  **[...] = momaperture_allring(..., 'ntimessplit',nsplit)**
   |          The number of recursive calls reducing the step size. Default 2
   | 
   |  **[...] = momaperture_allring(..., 'splitdivisor',splittdivisor)**
   |          The step divisor every time we split the step. Default 10.
   | 
   |  **[...] = momaperture_allring(..., 'verbose',verbose)**
   |          Print info the current position. Default 0.
   |          If set to 1 it will print info at every reference point.
   |          If set to 2 it will print info at each energy step.
   | 
   |  ex: **[map_l,map_h] = momaperture_allring(thering,points,nturns)**;
   |  ex: [map_l,map_h] = ...
   |    **momaperture_allring(thering,points,nturns,'reforbit',findorbit6(thering,points))**;

.. py:function:: TLT_IntPiw

   | integral in Piwinski Formula for the Lifetime

.. py:function:: TLT_IntPiw_k

   | integral in Piwinski Formula for the Lifetime with u=tan^2(k)

.. py:function:: TouschekPiwinskiLifeTime(..., 'abstol', abstol)

   | function [Tl,contributionsTL]=TouschekPiwinskiLifeTime(ring,dpp,Ib,...)
   | 
   |  evaluates Touschek Lifetime using Piwinski formula
   | 
   |  INPUT
   | 
   |  **touschekpiwinskilifetime**(
   |   ring,
   |   momentumaperturecolumnvector,  column array (size of r or positions)
   |                                  it can be length(r)x2, for positive and
   |                                  negative aperture
   |   current per bunch in A,        scalar
   |   positions where to evaluate,	  default: all elements with length>0  column array
   |   emittancex,                    default: atx modemittance(1)   scalar
   |   emittancey,                    default: emittancex/2		     scalar
   |   integration_method,            default: 'integral', may be: 'integral', 'quad', 'trapz')
   |   energy_spread,                 scalar
   |   bunch_length,	              scalar
   | 			   )
   | 
   |  OUTPUT
   | 
   |   contributionsTL 1/T contribution at each element
   |                   if dpp has positive and negative apertures, then contributionTL is a 2 columns vector
   | 
   |   Tl  Lifetime in seconds 1/Tl=sum(contributionsTL.*L)/sum(L);
   | 
   |  NAME-VALUE PAIRS
   | 
   |  **touschekpiwinskilifetime(..., 'abstol', abstol)**
   |    Absolute tolerance for the 'integral' function. Default: 1.0e-16
   | 
   |  **touschekpiwinskilifetime(..., 'reltol', abstol)**
   |    Relative tolerance for the 'integral' function. Default: 1.0e-12
   | 
   |  "The Touscheck Effect in strong focusing storage rings"
   |  A.Piwinski, DESY 98-179, November 1998
   | 
   |  "Touscheck Effect calculation and its applications to a transport line"
   |  A.Xiao M. Borland, Proceedings of PAC07, Albuquerque, New Mexico, USA
   | 

.. py:function:: momentum_aperture_at


