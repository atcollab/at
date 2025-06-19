CollectiveEffects
=================

.. rubric:: Functions


atbeam ATBEAM generates a particle distribution according to a sigma matrix
atsigma ATSIGMA constructs a beam sigma matrix 2x2 4x4 or 6x6

.. py:function:: atbeam

   |   PARTICLES=ATBEAM(NP,SIGMA)  Generate a particle distribution according to a sigma matrix
   |   PARTICLES=ATBEAM(NP,SIGMA,ORBIT) adds a center of mass to the distribution%
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

.. py:function:: atsigma

   | 
   |    SIGMA=ATSIGMA(BETA,ALPHA,EMIT)
   |        builds a 2x2 sigma matrix for a transverse plane
   | 
   |    SIGMA=ATSIGMA(ESPREAD,BLENGTH)
   |        builds a 2x2 sigma matrix for the longitudinal plane
   | 
   |    SIGMA=ATSIGMA(BETAX,ALPHAX,EMITX,BETAZ,ALPHAZ,EMITZ)
   |        builds a 4x4 sigma matrix
   | 
   |    SIGMA=ATSIGMA(BETAX,ALPHAX,EMITX,BETAZ,ALPHAZ,EMITZ,ESPREAD,BLENGTH)
   |        builds a 6x6 sigma matrix
   | 
   |    SIGMA=ATSIGMA(ATSTRUCT)
   |        builds a 6x6 sigma matrix
   | 
   |   See also atx

