.. _atphysics_module:

atphysics
=========

.. py:module:: atphysics

   Accelerator physics

.. toctree::
   :hidden:

   atphysics.CollectiveEffects
   atphysics.LinearOptics
   atphysics.LongitudinalDynamics
   atphysics.NonLinearDynamics
   atphysics.Orbit
   atphysics.ParameterSummaryFunctions
   atphysics.Radiation
   atphysics.TouschekPiwinski
   atphysics.TuneAndChromaticity
   atphysics.nafflib

.. rubric:: Modules


.. list-table::

   * - :ref:`collectiveeffects_module`
     - Collective Effects
   * - :ref:`linearoptics_module`
     - Linear beam dynamics
   * - :ref:`longitudinaldynamics_module`
     - Longitudinal beam dynamics
   * - :ref:`nonlineardynamics_module`
     - Non-linear beam dynamics
   * - :ref:`orbit_module`
     - Closed orbit search
   * - :ref:`parametersummaryfunctions_module`
     - Summaries of lattice parameters
   * - :ref:`radiation_module`
     - Functions dealing with synchrotron radiation
   * - :ref:`touschekpiwinski_module`
     - Touschek liftime and apertures
   * - :ref:`tuneandchromaticity_module`
     - Tune and chromaticity
   * - :ref:`nafflib_module`
     - NAFF frequency analysis

.. rubric:: Classes


.. list-table::

   * - :class:`PhysConstant`
     - Physical constants

.. rubric:: Functions


.. list-table::

   * - :func:`atavedata`
     - Average of optical functions on selected elements
   * - :func:`findrespm`
     - computes the change in the closed orbit due to parameter perturbations
   * - :func:`findspos`
     - returns longitudinal positions of accelerator lattice elements.

.. py:class:: PhysConstant

   | Physical constants
   
   | Automatically generated from https://physics.nist.gov/cuu/Constants/Table/allascii.txt
   
   | The values of the constants provided at this site are recommended for international use
   | by CODATA and are the latest available.
   | Termed the "2022 CODATA recommended values," they are generally recognized worldwide
   | for use in all fields of science and technology.
   | The values became available on 20 May 2024 and replaced the 2018 CODATA set.
   | They are based on all of the data available through 31 December 2022.
   | The 2022 adjustment was carried out under the auspices of the CODATA Task Group
   | on Fundamental Constants.

.. py:function:: atavedata(ring,dpp,refpts)

   | Average of optical functions on selected elements
   
   | **[lindata,avebeta,avemu,avedisp,tunes,chroms]=atavedata(ring,dpp,refpts)**
   
   | LINDATA : Identical to ATLINOPT output
   | AVEBEA :  Average Beta functions
   | AVEMU :   Average phase advance
   | AVEDISP : Average dispersion
   | TUNES : Vector of tunes
   | CHROMS : Vector of chromaticites
   
   | **[lindata,avebeta,avemu,avedisp,tunes,chroms]=atavedata(ring,dpp,refpts,orbitin)**
   |     does not search for closed orbit. instead ORBITIN is used

.. py:function:: findrespm(ring, obsindex, perturbindex, pvalue, 'field', m, n, orbitfunction, args)

   | computes the change in the closed orbit due to parameter perturbations
   |  Two calling syntax options
   |  1. **findrespm(ring, obsindex, perturbindex, pvalue, 'field', m, n, orbitfunction, args)**
   |  2. !!! not implemented yet **findrespm(ring, obsindex, perturbgroup, pvalue, orbitfunction, args)**
   
   |  RING      - ring lattice
   |  OBSINDEX  - indexes of elements where the orbit is observed (at the entrance)
   |  PERTURBINDEX  - Integer indexes of elements whose parameters are perturbed
   |                  used with syntax 1 only.
   
   |  PERTURBGROUP  - cell array of AT paramgroups. See ATPARAMGROUP
   |                used with syntax 2 only
   
   |  PVALUE    - amount of peturbation
   |              (Numeric array or scalar if all perturbations are the same magnitude)
   
   |  FIELD,M,N are only use with syntax 1.
   
   |  FIELD     - field name of the parameter to perturb (string)
   
   |  M,N       - index in the matrix, if the field is a matrix
   |              For example to perturb the quadrupole field in a
   |              multipole element
   |              FIELD = 'PolynomB', M = 1, N = 2
   
   |  ORBITFUNCTION  - specifies which of the FINDORBIT functions is used
   
   |              'findorbit4' (default)
   |              'findsyncorbit'
   |              'findorbit6'
   
   |  ARGS - additioanl arguments may be passsed to some of the FINDORBIT functions
   |              findorbit4     - constant momentum error dP
   |              findsyncorbit  - fixed orbit lengthening dCT
   
   
   |  Returns a 1-by-4 cell array of O-by-P matrixes
   |  where O = length(OBSINDEX) and P = length(PERTURB)
   |  one for each of the close orbit components: X, PX, Y, PY
   | See also :func:`atparamgroup`, :func:`findorbit`, :func:`findorbit4`, :func:`findorbit6`, :func:`findsyncorbit`

.. py:function:: findspos(ring,1:length(ring))

   | returns longitudinal positions of accelerator lattice elements.
   |   Return value is a row vector of positions S at the entrance of each
   |   element specified REFPTS (index list or logical mask)
   
   |  Note: REFPTS is an array of increasing indexes that
   |    select elements from range 1 to length(LATTICE)+1.
   |    REFPTS is allowed to go 1 point beyond the
   |    number of elements. In this case the last point is
   |    the EXIT of the last element. If LATTICE is a RING
   |    it is also the entrance of the first element after 1 turn.
   
   |    Note:
   |    1. Use **findspos(ring,1:length(ring))** for to find
   |       longitudinal position of all elements
   |    2. Use **findspos(line,length(line)+1)** to find the
   |       total physical length
   |    3. If line is a closed ring, exit of the last element
   |       is also the entrance to the first.

