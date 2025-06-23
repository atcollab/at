.. _collectiveeffects_module:

CollectiveEffects
=================

.. rubric:: Functions


.. list-table::

   * - :func:`atbeam`
     - generates a particle distribution according to a sigma matrix
   * - :func:`atsigma`
     - constructs a beam sigma matrix 2x2 4x4 or 6x6

.. py:function:: atbeam(np,sigma)

   | generates a particle distribution according to a sigma matrix
   |   **particles=atbeam(np,sigma)**  Generate a particle distribution according to a sigma matrix
   |   **particles=atbeam(np,sigma,orbit)** adds a center of mass to the distribution%
   | 
   |   INPUTS
   |     1. NP     number of particles
   |     2. SIGMA  beam matrix (2x2, 4x4, 6x6)
   |     3. ORBIT  closed orbit
   | 
   |   OUPUTS
   |     1. PARTICLES particle distribution
   | 
   |   NOTES
   |     1. random generator is randn
   | 
   |   See also atsigma

.. py:function:: atsigma(beta,alpha,emit)

   | constructs a beam sigma matrix 2x2 4x4 or 6x6
   | 
   |    **sigma=atsigma(beta,alpha,emit)**
   |        builds a 2x2 sigma matrix for a transverse plane
   | 
   |    **sigma=atsigma(espread,blength)**
   |        builds a 2x2 sigma matrix for the longitudinal plane
   | 
   |    **sigma=atsigma(betax,alphax,emitx,betaz,alphaz,emitz)**
   |        builds a 4x4 sigma matrix
   | 
   |    **sigma=atsigma(betax,alphax,emitx,betaz,alphaz,emitz,espread,blength)**
   |        builds a 6x6 sigma matrix
   | 
   |    **sigma=atsigma(atstruct)**
   |        builds a 6x6 sigma matrix
   | 
   |   See also atx

