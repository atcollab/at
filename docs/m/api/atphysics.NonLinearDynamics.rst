.. _nonlineardynamics_module:

NonLinearDynamics
=================

.. rubric:: Functions


.. list-table::

   * - :func:`RDTbuildupFluct`
     - 
   * - :func:`RDTfluctuationIndicator`
     - quantitative representation of RDT fluctuations
   * - :func:`atnuampl`
     - computes tune shift with amplitude
   * - :func:`computeRDT`
     - Computes Hamiltonian resonance driving terms (RDTs)
   * - :func:`computeRDTfluctuation`
     - Computes Hamiltonian resonance driving terms (RDTs)
   * - :func:`tunespaceplot`
     - draws a tune diagram

.. py:function:: RDTbuildupFluct


.. py:function:: RDTfluctuationIndicator(ring,varargin)

   | quantitative representation of RDT fluctuations
   |    This function calls computeRDTfluctuation(ring, varargin) to compute
   |    RDT fluctuations, and provides one example to quantitatively represents
   |    the RDT fluctuations.
   
   |  **[h3ave,h4ave,h3chroave,adts]=RDTfluctuationIndicator(ring,varargin)**
   
   |    ring is the AT lattice
   |    The additional argument:
   |      nslices: number of slices of each sextupole, which affects the
   |               crossing terms. default: 4.
   
   |    h3ave and h4ave quantitatively represents 3rd- and 4th-order geometric
   |      RDT fluctuations, respectively.
   
   |    h3ave + w * h4ave can be used to represents geometric RDT fluctuations.
   |      The coefficient w can be estimated by action of particle at the
   |      anticipated DA, w ~ (2J_x)^0.5 = x / betax^0.5,  usually 0.01 is OK.
   
   |    h3chroave is the fluctuation of 3rd-order chromatic RDTs,
   |      defined similarly to h3ave.
   
   |    adts is the sum of three ADTS terms, which are also used in onlinear optimization.
   |    The fluctuations of ADTS terms are not considered.
   |    It is calculated here to avoid duplicate computation.
   
   |  Noting:
   |     1.This function provides one example to quantitatively represents the
   |       RDT fluctuations similar to that in Ref.[1]. But in Ref.[1], only
   |       the RDTs at the locations of nonlinear magnets are considered.
   |       We think the differences are small and it's more important to
   |       keep the function simple.
   |     2.People can call computeRDTfluctuation(ring, varargin) directly,
   |       and try other quantitative representations.
   |     3.The build-up RDT fluctuation can also be used, see Ref.[1].
   |       If the build-up RDT fluctuations are used, it is better to calculate
   |       the RDT build-up fluctuations for multiple periods to have
   |       better convergence of calculation.
   
   |  REFERENCE:
   |    [1] B. Wei, Z. Bai, J. Tan, L. Wang, and G. Feng, Phys. Rev. Accel. Beams 26, 084001 (2023)
   

.. py:function:: atnuampl(ring,amplitude)

   | computes tune shift with amplitude
   | **[nux,nuz]=atnuampl(ring,amplitude)**
   | **[nux,nuz]=atnuampl(ring,amplitude,1)**
   
   | 	Computes tunes for the specified horizontal amplitudes
   
   | **[nux,nuz]=atnuampl(ring,amplitude,3)**
   
   | 	Computes tunes for the specified vertical amplitudes
   
   | **atnuampl(...)**
   |    Plots the computed tunes in the current axes
   
   | **atnuampl(...,name,value)**
   |    Uses additional options specified by one or more Name,Value pairs.
   |    Possible values are:
   |        orbit:  initial closed orbit
   |        nturns: specify the number of turns for tracking (default 256)
   |        method: specify the method for tune determination
   |                1: Highest peak in fft
   |                2: Interpolation on fft results
   |                3: Windowing + interpolation (default)
   |                4: NAFF
   |    Other options are transmitted to the plot function

.. py:function:: computeRDT(ring, index, varargin)

   | Computes Hamiltonian resonance driving terms (RDTs)
   |    This function calls RDTElegantAT mex function and returns the
   |    hamiltonian resonance driving terms, using the elegant c++
   |    function computeDrivingTerms().
   
   |    **rdt=computeRDT(ring, index, varargin)**
   
   |    ring is the AT lattice
   |    index is the vector of indexes where one wants to compute RDTs
   |    The additional arguments can be up to five strings:
   |    chromatic, coupling, geometric1, geometric2 and tuneshifts
   
   |    example:
   |    **rdt=computeRDT(ring, indexbpm, 'geometric1', 'tuneshifts')**;
   |    creates an array of structs (the length of the array is the number of
   |    indexes where you want to compute driving terms) with first order
   |    geometric driving terms and tune shifts with amplitude.
   |    The driving terms are complex numbers, the tune shifts are real.
   

.. py:function:: computeRDTfluctuation(ring, varargin)

   | Computes Hamiltonian resonance driving terms (RDTs)
   |    This function is based on simplestoragering and returns the RDTs
   |    and their longitudinal fluctuations.
   
   |  **[rdt,buildup_fluctuation,natural_fluctuation]=computeRDTfluctuation(ring, varargin)**
   
   |    ring is the AT lattice
   |   The additional arguments:
   |    nslices: number of slices of each sextupole, which affects the crossing
   |        terms. default: 4.
   |    nperiods: number of periods. RDTs and RDT build-up fluctuations will be
   |        computed for n periods.  default: 1.
   |        natural RDT fluctuation of different periods are the same.
   |        So the results contain only one period.
   
   |    RDT: struct, RDTs (complex numbers) and
   |        amplitude-dependent tune shifts (real)
   |        (ADTS are calculated using h22000, h11110 and h00220)
   |    buildup_fluctuation: a struct of complex arrays,
   |        accumulated RDTs from s=0,
   |        showing the build-up and cancellation of RDTs along the
   |        longitudinal position.
   |    natural_fluctuation: a struct of double arrays,
   |        absolute values of one-period RDTs observed at different
   |        longitudinal starting position.
   |        same as s_dependent_driving_terms in ELEGANT.
   
   |   REFERENCES
   |     [1] Johan Bengtsson, SLS Note 09/97, (1997)
   |     [2] S. C. Leemann, A. Streun, Phys. Rev. ST Accel. Beams 14, 030701 (2011)
   |     [3] A. Franchi, L. Farvacque, F. Ewald, G. Le Bec, and K. B. Scheidt, Phys. Rev. ST Accel. Beams 17, 074001 (2014)
   |     [4] B. Wei, Z. Bai, J. Tan, L. Wang, and G. Feng, Phys. Rev. Accel. Beams 26, 084001 (2023)
   

.. py:function:: tunespaceplot

   | draws a tune diagram
   |  resonance lines: m*nu_x + n*nu_y = p

