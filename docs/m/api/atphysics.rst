atphysics
=========

.. toctree::
   :hidden:

   Orbit <atphysics.Orbit>
   nafflib <atphysics.nafflib>
   Radiation <atphysics.Radiation>
   LinearOptics <atphysics.LinearOptics>
   NonLinearDynamics <atphysics.NonLinearDynamics>
   TuneAndChromaticity <atphysics.TuneAndChromaticity>
   LongitudinalDynamics <atphysics.LongitudinalDynamics>
   ParameterSummaryFunctions <atphysics.ParameterSummaryFunctions>
   TouschekPiwinski <atphysics.TouschekPiwinski>
   CollectiveEffects <atphysics.CollectiveEffects>

.. rubric:: Modules


Orbit ORBIT
nafflib NAFFLIB
Radiation RADIATION
LinearOptics LINEAROPTICS
NonLinearDynamics NONLINEARDYNAMICS
TuneAndChromaticity TUNEANDCHROMATICITY
LongitudinalDynamics LONGITUDINALDYNAMICS
ParameterSummaryFunctions PARAMETERSUMMARYFUNCTIONS
TouschekPiwinski TOUSCHEKPIWINSKI
CollectiveEffects COLLECTIVEEFFECTS

.. rubric:: Functions


findspos FINDSPOS returns longitudinal positions of accelerator lattice elements.
atavedata ATAVEDATA       Average of optical functions on selected elements
PhysConstant Physical constants
findrespm FINDRESPM computes the change in the closed orbit due to parameter perturbations

.. py:function:: findspos

   |   Return value is a row vector of positions S at the entrance of each
   |   element specified REFPTS (index list or logical mask)
   | 
   |  Note: REFPTS is an array of increasing indexes that
   |    select elements from range 1 to length(LATTICE)+1.
   |    REFPTS is allowed to go 1 point beyond the
   |    number of elements. In this case the last point is
   |    the EXIT of the last element. If LATTICE is a RING
   |    it is also the entrance of the first element after 1 turn.
   | 
   |    Note:
   |    1. Use findspos(RING,1:length(RING)) for to find
   |       longitudinal position of all elements
   |    2. Use findspos(LINE,length(LINE)+1) to find the
   |       total physical length
   |    3. If line is a closed ring, exit of the last element
   |       is also the entrance to the first.

.. py:function:: atavedata

   | 
   | [LINDATA,AVEBETA,AVEMU,AVEDISP,TUNES,CHROMS]=ATAVEDATA(RING,DPP,REFPTS)
   | 
   | LINDATA : Identical to ATLINOPT output
   | AVEBEA :  Average Beta functions
   | AVEMU :   Average phase advance
   | AVEDISP : Average dispersion
   | TUNES : Vector of tunes
   | CHROMS : Vector of chromaticites
   | 
   | [LINDATA,AVEBETA,AVEMU,AVEDISP,TUNES,CHROMS]=ATAVEDATA(RING,DPP,REFPTS,ORBITIN)
   |     does not search for closed orbit. instead ORBITIN is used

.. py:function:: PhysConstant

   | 
   | Automatically generated from https://physics.nist.gov/cuu/Constants/Table/allascii.txt
   | 
   | The values of the constants provided at this site are recommended for international use
   | by CODATA and are the latest available.
   | Termed the "2022 CODATA recommended values," they are generally recognized worldwide
   | for use in all fields of science and technology.
   | The values became available on 20 May 2024 and replaced the 2018 CODATA set.
   | They are based on all of the data available through 31 December 2022.
   | The 2022 adjustment was carried out under the auspices of the CODATA Task Group
   | on Fundamental Constants.

.. py:function:: findrespm

   |  Two calling syntax options 
   |  1. FINDRESPM(RING, OBSINDEX, PERTURBINDEX, PVALUE, 'FIELD', M, N, ORBITFUNCTION, ARGS)
   |  2. !!! not implemented yet FINDRESPM(RING, OBSINDEX, PERTURBGROUP, PVALUE, ORBITFUNCTION, ARGS)
   | 
   |  RING      - ring lattice
   |  OBSINDEX  - indexes of elements where the orbit is observed (at the entrance)
   |  PERTURBINDEX  - Integer indexes of elements whose parameters are perturbed
   |                  used with syntax 1 only. 
   |              
   |  PERTURBGROUP  - cell array of AT paramgroups. See ATPARAMGROUP
   |                used with syntax 2 only
   | 
   |  PVALUE    - amount of peturbation 
   |              (Numeric array or scalar if all perturbations are the same magnitude) 
   |  
   |  FIELD,M,N are only use with syntax 1. 
   | 
   |  FIELD     - field name of the parameter to perturb (string)
   | 
   |  M,N       - index in the matrix, if the field is a matrix
   |              For example to perturb the quadrupole field in a
   |              multipole element
   |              FIELD = 'PolynomB', M = 1, N = 2
   | 
   |  ORBITFUNCTION  - specifies which of the FINDORBIT functions is used
   |              
   |              'findorbit4' (default)
   |              'findsyncorbit'
   |              'findorbit6'
   |              
   |  ARGS - additioanl arguments may be passsed to some of the FINDORBIT functions
   |              findorbit4     - constant momentum error dP
   |              findsyncorbit  - fixed orbit lengthening dCT
   |               
   | 
   |  Returns a 1-by-4 cell array of O-by-P matrixes 
   |  where O = length(OBSINDEX) and P = length(PERTURB)
   |  one for each of the close orbit components: X, PX, Y, PY
   |  See also ATPARAMGROUP, FINDORBIT, FINDORBIT4, FINDORBIT6, FINDSYNCORBIT

