.. _plotfunctions_module:

plotfunctions
=============

.. rubric:: Functions


.. list-table::

   * - :func:`CurlyH`
     - function [H,Hv]=CurlyH(RING,dp,ind)
   * - :func:`CurlyHlindata`
     - function [H,Hv]=CurlyHlindata(lindata)
   * - :func:`plBeamSize`
     - Plot H and V beam size
   * - :func:`plClosedOrbit`
     - Plots H and V 4x4 closed orbit
   * - :func:`plCorrectorStrength`
     - Plot PolynomB
   * - :func:`plEmitContrib`
     - Plot H/rho³ at every dipole
   * - :func:`plPolynomBComp`
     - PlotBn coefficient with normalization
   * - :func:`plPolynomBSxtOct`
     - Plots Bn for sextupole and octupole magnets
   * - :func:`plSigmaSigmap`
     - Plots beam sizes and divergences
   * - :func:`plenvelope`
     - Plot beam envelope
   * - :func:`plotAperture`
     - Plots x and y aperture.
   * - :func:`plotB0curlyh`
     - Plot B and H
   * - :func:`plotERAperture`
     - Plot RApertures EApertures
   * - :func:`plotRDT`
     - 
   * - :func:`plotWdispP`
     - plotWdispP    Plot W functions
   * - :func:`plot_betabeat`
     - function plot_betabeat(THERING_ref,THERING_mod)
   * - :func:`plot_trajectory`
     - Plots particle trajectories
   * - :func:`plotbetadisp`
     - function [s,plotdata]=plotbetadisp(ring,dpp,plotfun,varargin)
   * - :func:`plotbetadispcurlyh`
     - Plot beta, dispersion and H
   * - :func:`plotsqrtbetadispcurlyh`
     - Plot sqrt(beta), dispersion and H
   * - :func:`pltouschek`
     - Plots Touschek lifetime contribution
   * - :func:`plxi`
     - plotdata=plxi(lindata,ring,dpp)

.. py:function:: CurlyH

   | function [H,Hv]=CurlyH(RING,dp,ind)
   
   |  computes Curly H (dispersion invariant)
   
   |  RING :at lattice
   |  dp  : energy deviation
   |  ind : reference positions
   
   |  output
   
   |  H  : horizontal dispersion invariant
   |  Hv : vertical dispersion invariant
   
   

.. py:function:: CurlyHlindata

   | function [H,Hv]=CurlyHlindata(lindata)
   
   |  computes Curly H (dispersion invariant)
   
   |  lindata is the first ouptut of [lindata,~,~]=atlinopt(...)
   |  (include at least 2 ouptut arguments for dispersion computation)
   
   |  output
   
   |  H  : horizontal dispersion invariant
   |  Hv : vertical dispersion invariant
   
   

.. py:function:: plBeamSize

   | Plot H and V beam size
   
   |  USAGE:
   |  >> atbaseplot(ring,@**plBeamSize**);
   |  >> atplot(ring,@**plBeamSize**);     (obsolete)
   
   | See also :func:`atbaseplot`

.. py:function:: plClosedOrbit

   | Plots H and V 4x4 closed orbit
   
   | Helper function for atplot: plot
   | - H and V closed orbits on left axis
   | - Dispersion on right axis
   
   |   EXAMPLEs
   |  >> atbaseplot(ring,@**plClosedOrbit**,{'dp',0.01});
   |  >> atplot(ring,@**plClosedOrbit**,'dp',0.01);     (obsolete)
   
   | See also :func:`atplot`, :func:`atbaseplot`

.. py:function:: plCorrectorStrength

   | Plot PolynomB
   | Helper function for atplot: plot
   | - PolynomB(1), PolynomA(1), PolynomA(2), PolynomA(2) on the left axis
   
   |   EXAMPLEs
   |  >> atbaseplot(ring,'synopt',false,@**plCorrectorStrength**);
   |  >> atplot(ring,@**plCorrectorStrength**,'synopt',false);     (obsolete)
   
   
   | See also :func:`atplot`, :func:`atbaseplot`

.. py:function:: plEmitContrib

   | Plot H/rho³ at every dipole
   
   |  USAGE:
   |  >> atbaseplot(ring,@**plEmitContrib**);
   |  >> atplot(ring,@**plEmitContrib**);     (obsolete)
   
   | See also :func:`atbaseplot`

.. py:function:: plPolynomBComp

   | PlotBn coefficient with normalization
   |  **plPolynomBComp** default plotting function for ATPLOT
   
   |  Plots polynomB for ring and ring1

.. py:function:: plPolynomBSxtOct

   | Plots Bn for sextupole and octupole magnets
   | DEFAULTPLOT    Default plotting function for ATPLOT
   
   | Plots polynomB for ring and ring1

.. py:function:: plSigmaSigmap

   | Plots beam sizes and divergences
   |  Plots sigmax and sigmay on left axis and
   |        sigmax' and sigmay' on right axis

.. py:function:: plenvelope

   | Plot beam envelope
   
   | Helper function for atplot: plot
   | - H and V beam envelopes on left axis
   
   |  USAGE:
   |  >> atbaseplot(ring,@**plenvelope**);
   |  >> atplot(ring,@**plenvelope**);     (obsolete)
   
   | See also :func:`atbaseplot`

.. py:function:: plotAperture

   | Plots x and y aperture.
   
   | Helper function for atplot: plot the physical aperture
   
   |   USAGE:
   |  >> atbaseplot(ring,@**plotAperture**);
   |  >> atplot(ring,@**plotAperture**);        (obsolete)
   
   | See also :func:`atplot`, :func:`atbaseplot`

.. py:function:: plotB0curlyh

   | Plot B and H
   
   |  USAGE:
   |  >> atbaseplot(ring,@**plotB0curlyh**);
   |  >> atplot(ring,@**plotB0curlyh**);     (obsolete)
   
   | See also :func:`atbaseplot`

.. py:function:: plotERAperture

   | Plot RApertures EApertures
   
   | Helper function for atplot:
   |  plot the Elliptic and Rectangular physical apertures
   
   |   USAGE:
   |  >> atbaseplot(ring,@**plotERAperture**);
   |  >> atplot(ring,@**plotERAperture**);      (obsolete)
   
   | See also :func:`atplot`, :func:`atbaseplot`

.. py:function:: plotRDT

   
   |   **plotRDT** plots the absolute value of the hamiltonian terms
   |   **plotRDT** must be used with atplot:
   
   |   atplot(ring,@**plotRDT**,'geometric1') plots the first order geometric terms
   |   atplot(ring,@**plotRDT**,'chromatic') plots the chromatic terms
   |   atplot(ring,@**plotRDT**,'geometric2') plots the second order geometric terms
   |   atplot(ring,@**plotRDT**,'coupling') plots the coupling terms
   
   |   see also: computeRDT, atplot

.. py:function:: plotWdispP

   | plotWdispP    Plot W functions
   
   | Helper function for atplot: plot
   | - W functions (derivatives of beta-functions versus momentum) on left axis
   | - derivative of dispersion on right axis

.. py:function:: plot_betabeat

   | function plot_betabeat(THERING_ref,THERING_mod)
   
   |  returns plot of beta beat of THERING_mod respect to THERING_ref

.. py:function:: plot_trajectory

   | Plots particle trajectories
   
   | Helper function for atplot: plot
   | - H and V trajectories on the left axis
   
   |  USAGE:
   |  >> atbaseplot(ring,@**plot_trajectory**,{INPUT_COORDS);
   |  >> atplot(ring,@**plot_trajectory**,INPUT_COORDS);     (obsolete)
   
   | See also :func:`atbaseplot`

.. py:function:: plotbetadisp

   | function [s,plotdata]=plotbetadisp(ring,dpp,plotfun,varargin)
   | **plotbetadisp** Plot beta functions and dispersion
   
   | Helper function for tplot:
   | - beta functions on the left axis
   | - diespersion on the right axis
   
   |   USAGE:
   |  >> atbaseplot(ring,@**plotbetadisp**);
   |  >> atplot(ring,@**plotbetadisp**);      (obsolete)
   
   | See also :func:`atplot`, :func:`atbaseplot`

.. py:function:: plotbetadispcurlyh

   | Plot beta, dispersion and H
   
   |  USAGE:
   |  >> atbaseplot(ring,@**plotbetadispcurlyh**);
   |  >> atplot(ring,@**plotbetadispcurlyh**);     (obsolete)
   
   | See also :func:`atbaseplot`

.. py:function:: plotsqrtbetadispcurlyh

   | Plot sqrt(beta), dispersion and H
   
   |  USAGE:
   |  >> atbaseplot(ring,@**plotsqrtbetadispcurlyh**);
   |  >> atplot(ring,@**plotsqrtbetadispcurlyh**);     (obsolete)
   
   | See also :func:`atbaseplot`

.. py:function:: pltouschek(lindata,ring,dpp)

   | Plots Touschek lifetime contribution
   | **plotdata=pltouschek(lindata,ring,dpp)**
   |  plots curly H function and also 1/(sigx sigy) to understand local scattering rate

.. py:function:: plxi

   | plotdata=plxi(lindata,ring,dpp)
   | xi function in Touschek formula gives trans. velocity in beam frame

