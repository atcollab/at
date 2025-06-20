.. _tuneandchromaticity_module:

TuneAndChromaticity
===================

.. rubric:: Functions


.. list-table::

   * - :func:`fitchrom2`
     - FITCHROM2 fits chromaticity  of THERING using 2 sextupole families
   * - :func:`findtune`
     - FINDTUNE   get the tune value from turn by turn positions
   * - :func:`fittune2`
     - FITTUNE2 fits linear tunes of THERING using 2 quadrupole families
   * - :func:`tunechrom`
     - TUNECHROM computes linear tunes and chromaticities
   * - :func:`intfft`
     - INTFFT Calculates the tune from interpolated FFT of the trajectory.

.. py:function:: fitchrom2

   |  FITCHROM2(NEWCHROM,SEXTUPOLEFAMILY1,SEXTUPOLEFAMILY2)

.. py:function:: findtune

   | 
   | TUNE=FINDTUNE(POS,METHOD)
   | 
   | POS:       Tune-by-turn particle position
   | METHOD:    Method for tune determination:
   |                1: Highest peak in fft
   |                2: Interpolation on fft results
   |                3: Windowing + interpolation
   | 
   | [TUNE,SPECTRUM]=FINDTUNE(...) Also returns the fft

.. py:function:: fittune2

   |  FITTUNE2(NEWTUNES,QUADFAMILY1,QUADFAMILY2)
   |  [dK, Jinv] = FITTUNE2(NEWTUNES,QUADFAMILY1,QUADFAMILY2)

.. py:function:: tunechrom

   | 
   | TUNE=TUNECHROM(RING)	Quick calculation of the fractional part of the tune
   | 	from numerically computed transfer matrix.
   | 
   |  [TUNE, CHROM] = TUNECHROM(RINGD,DP,'get_chrom') - optionally computes the
   |     chromaticities by numerical differentiation from the difference between
   |    tune values at momentums DP+0.5*DPStep and DP-0.5*DPStep
   | 
   | [...]=TUNECHROM(...,'orbit',ORBITIN	Do not search for closed orbit.
   |    Instead ORBITIN,a 6x1 vector of initial conditions is used:
   |    This syntax is useful to avoid recomputing the closed orbit if is
   |    already known;
   | 
   | [...]=TUNECHROM(RING,DP)       (obsolete)
   | [...]=TUNECHROM(RING,...,'dp',DP)	Specify the momentum deviation when
   |    radiation is OFF (default: 0)
   | 
   | [...]=TUNECHROM(RING,...,'dct',DCT) Specify the path lengthening when
   |    radiation is OFF (default: 0)
   | 
   | [...]=TUNECHROM(RING,...,'df',DF) Specify the RF frequency deviation when
   |    radiation is OFF (default: 0)
   | 
   |  Note: TUNECHROM computes tunes and chromaticities from the one-turn
   |    transfer matrix. The transfer matrix is computed from tracking using
   |    numerical differentiation. The error of numerical differentiation
   |    is sensitive to the step size. (Reference: Numerical Recipes)
   |    The calculation of tunes involves one numerical differentiation.
   |    The calculation of chromaticity involves TWO!!! numerical differentiations.
   |    The error in calculated chromaticity from may be substantial (~ 1e-5).
   |    Use the XYStep and DPStep keyword arguments to control the step size
   |    in chromaticity calculations
   | 
   |  See also ATLINOPT6

.. py:function:: intfft

   |  INTFFT(X) X must be a column vector.
   |   If X is a matrix - each column is treated as
   |   a separate trajectory
   |  INTFFT(X,GUESS,DELTA) searches for peaks in the FFT spectrum
   |   only within the range (X-DELTA ... X+DELTA) 
   |   The same values of GUESS and DELTA are used for all columns of X 

