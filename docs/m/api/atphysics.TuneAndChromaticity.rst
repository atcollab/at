.. _tuneandchromaticity_module:

TuneAndChromaticity
===================

.. rubric:: Functions


.. list-table::

   * - :func:`findtune`
     - get the tune value from turn by turn positions
   * - :func:`fitchrom2`
     - fits chromaticity  of THERING using 2 sextupole families
   * - :func:`fittune2`
     - fits linear tunes of THERING using 2 quadrupole families
   * - :func:`intfft`
     - Calculates the tune from interpolated FFT of the trajectory.
   * - :func:`tunechrom`
     - computes linear tunes and chromaticities

.. py:function:: findtune(pos,method)

   | get the tune value from turn by turn positions
   | 
   | **tune=findtune(pos,method)**
   | 
   | POS:       Tune-by-turn particle position
   | METHOD:    Method for tune determination:
   |                1: Highest peak in fft
   |                2: Interpolation on fft results
   |                3: Windowing + interpolation
   | 
   | **[tune,spectrum]=findtune(...)** Also returns the fft

.. py:function:: fitchrom2(newchrom,sextupolefamily1,sextupolefamily2)

   | fits chromaticity  of THERING using 2 sextupole families
   |  **fitchrom2(newchrom,sextupolefamily1,sextupolefamily2)**

.. py:function:: fittune2(newtunes,quadfamily1,quadfamily2)

   | fits linear tunes of THERING using 2 quadrupole families
   |  **fittune2(newtunes,quadfamily1,quadfamily2)**
   |  **[dk, jinv] = fittune2(newtunes,quadfamily1,quadfamily2)**

.. py:function:: intfft(x)

   | Calculates the tune from interpolated FFT of the trajectory.
   |  **intfft(x)** X must be a column vector.
   |   If X is a matrix - each column is treated as
   |   a separate trajectory
   |  **intfft(x,guess,delta)** searches for peaks in the FFT spectrum
   |   only within the range (X-DELTA ... X+DELTA)
   |   The same values of GUESS and DELTA are used for all columns of X

.. py:function:: tunechrom(ring)

   | computes linear tunes and chromaticities
   | 
   | **tune=tunechrom(ring)**	Quick calculation of the fractional part of the tune
   | 	from numerically computed transfer matrix.
   | 
   |  **[tune, chrom] = tunechrom(ringd,dp,'get_chrom')** - optionally computes the
   |     chromaticities by numerical differentiation from the difference between
   |    tune values at momentums DP+0.5*DPStep and DP-0.5*DPStep
   | 
   | **[...]=tunechrom**(...,'orbit',ORBITIN	Do not search for closed orbit.
   |    Instead ORBITIN,a 6x1 vector of initial conditions is used:
   |    This syntax is useful to avoid recomputing the closed orbit if is
   |    already known;
   | 
   | **[...]=tunechrom(ring,dp)       (obsolete)**
   | **[...]=tunechrom(ring,...,'dp',dp)**	Specify the momentum deviation when
   |    radiation is OFF (default: 0)
   | 
   | **[...]=tunechrom(ring,...,'dct',dct)** Specify the path lengthening when
   |    radiation is OFF (default: 0)
   | 
   | **[...]=tunechrom(ring,...,'df',df)** Specify the RF frequency deviation when
   |    radiation is OFF (default: 0)
   | 
   |  Note: **tunechrom** computes tunes and chromaticities from the one-turn
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

