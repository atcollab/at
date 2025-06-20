.. _plotfunctions_module:

plotfunctions
=============

.. rubric:: Functions


.. list-table::

   * - :func:`plotbetadispcurlyh`
     - PLOTBETADISPCURLYH Plot beta, dispersion and H
   * - :func:`plEmitContrib`
     - PLEMITCONTRIB  Plot H/rho³ at every dipole
   * - :func:`plxi`
     - plotdata=plxi(lindata,ring,dpp)
   * - :func:`plotAperture`
     - PLOTAPERTURE Plots x and y aperture.
   * - :func:`pltouschek`
     - PLTOUSCHEK Plots Touschek lifetime contribution
   * - :func:`plPolynomBSxtOct`
     - PLPOLYNOMBSXTOCT Plots Bn for sextupole and octupole magnets
   * - :func:`plSigmaSigmap`
     - PLSIGMASIGMAP Plots beam sizes and divergences
   * - :func:`plot_trajectory`
     - PLOT_TRAJECTORY    Plots particle trajectories
   * - :func:`plotERAperture`
     - PLOTERAPERTURE Plot RApertures EApertures
   * - :func:`plCorrectorStrength`
     - PLCORRECTORSTRENGTH Plot PolynomB
   * - :func:`plotsqrtbetadispcurlyh`
     - PLOTSQRTBETADISPCURLYH Plot sqrt(beta), dispersion and H
   * - :func:`plot_betabeat`
     - function plot_betabeat(THERING_ref,THERING_mod)
   * - :func:`plotRDT`
     - 
   * - :func:`plenvelope`
     - PLENVELOPE    Plot beam envelope
   * - :func:`CurlyHlindata`
     - function [H,Hv]=CurlyHlindata(lindata)
   * - :func:`plotbetadisp`
     - function [s,plotdata]=plotbetadisp(ring,dpp,plotfun,varargin)
   * - :func:`plPolynomBComp`
     - PLPOLYNOMBCOMP PlotBn coefficient with normalization
   * - :func:`plClosedOrbit`
     - PLCLOSEDORBIT Plots H and V 4x4 closed orbit
   * - :func:`plotB0curlyh`
     - PLOTB0CURLYH  Plot B and H
   * - :func:`CurlyH`
     - function [H,Hv]=CurlyH(RING,dp,ind)
   * - :func:`plBeamSize`
     - PLBEAMSIZE Plot H and V beam size
   * - :func:`plotWdispP`
     - plotWdispP    Plot W functions

.. py:function:: plotbetadispcurlyh

   | PLOTBETADISPCURLYH Plot beta, dispersion and H
   | 
   |  USAGE:
   |  >> atbaseplot(ring,@PLOTBETADISPCURLYH);
   |  >> atplot(ring,@PLOTBETADISPCURLYH);     (obsolete)
   | 
   |   See also atbaseplot

.. py:function:: plEmitContrib

   | PLEMITCONTRIB  Plot H/rho³ at every dipole
   | 
   |  USAGE:
   |  >> atbaseplot(ring,@PLEMITCONTRIB);
   |  >> atplot(ring,@PLEMITCONTRIB);     (obsolete)
   | 
   |   See also atbaseplot

.. py:function:: plxi

   | plotdata=plxi(lindata,ring,dpp)
   | xi function in Touschek formula gives trans. velocity in beam frame

.. py:function:: plotAperture

   | PLOTAPERTURE Plots x and y aperture.
   | 
   | Helper function for atplot: plot the physical aperture
   | 
   |   USAGE:
   |  >> atbaseplot(ring,@plotAperture);
   |  >> atplot(ring,@plotAperture);        (obsolete)
   | 
   | See also atplot atbaseplot

.. py:function:: pltouschek

   | PLTOUSCHEK Plots Touschek lifetime contribution
   | plotdata=pltouschek(lindata,ring,dpp)
   |  plots curly H function and also 1/(sigx sigy) to understand local scattering rate

.. py:function:: plPolynomBSxtOct

   | PLPOLYNOMBSXTOCT Plots Bn for sextupole and octupole magnets
   | DEFAULTPLOT    Default plotting function for ATPLOT
   | 
   | Plots polynomB for ring and ring1

.. py:function:: plSigmaSigmap

   | PLSIGMASIGMAP Plots beam sizes and divergences
   |  Plots sigmax and sigmay on left axis and
   |        sigmax' and sigmay' on right axis

.. py:function:: plot_trajectory

   | PLOT_TRAJECTORY    Plots particle trajectories
   | 
   | Helper function for atplot: plot
   | - H and V trajectories on the left axis
   | 
   |  USAGE:
   |  >> atbaseplot(ring,@PLOT_TRAJECTORY,{INPUT_COORDS);
   |  >> atplot(ring,@PLOT_TRAJECTORY,INPUT_COORDS);     (obsolete)
   | 
   |   See also atbaseplot

.. py:function:: plotERAperture

   | PLOTERAPERTURE Plot RApertures EApertures
   | 
   | Helper function for atplot:
   |  plot the Elliptic and Rectangular physical apertures
   | 
   |   USAGE:
   |  >> atbaseplot(ring,@plotERAperture);
   |  >> atplot(ring,@plotERAperture);      (obsolete)
   | 
   | See also atplot atbaseplot

.. py:function:: plCorrectorStrength

   | PLCORRECTORSTRENGTH Plot PolynomB
   | Helper function for atplot: plot
   | - PolynomB(1), PolynomA(1), PolynomA(2), PolynomA(2) on the left axis
   | 
   |   EXAMPLEs
   |  >> atbaseplot(ring,'synopt',false,@plCorrectorStrength);
   |  >> atplot(ring,@plCorrectorStrength,'synopt',false);     (obsolete)
   | 
   |   See also atplot atbaseplot
   | 

.. py:function:: plotsqrtbetadispcurlyh

   | PLOTSQRTBETADISPCURLYH Plot sqrt(beta), dispersion and H
   | 
   |  USAGE:
   |  >> atbaseplot(ring,@PLOTSQRTBETADISPCURLYH);
   |  >> atplot(ring,@PLOTSQRTBETADISPCURLYH);     (obsolete)
   | 
   |   See also atbaseplot

.. py:function:: plot_betabeat

   | function plot_betabeat(THERING_ref,THERING_mod)
   | 
   |  returns plot of beta beat of THERING_mod respect to THERING_ref

.. py:function:: plotRDT

   | 
   |   plotRDT plots the absolute value of the hamiltonian terms
   |   plotRDT must be used with atplot:
   | 
   |   atplot(ring,@plotRDT,'geometric1') plots the first order geometric terms
   |   atplot(ring,@plotRDT,'chromatic') plots the chromatic terms
   |   atplot(ring,@plotRDT,'geometric2') plots the second order geometric terms
   |   atplot(ring,@plotRDT,'coupling') plots the coupling terms
   | 
   |   see also: computeRDT, atplot

.. py:function:: plenvelope

   | PLENVELOPE    Plot beam envelope
   | 
   | Helper function for atplot: plot
   | - H and V beam envelopes on left axis
   | 
   |  USAGE:
   |  >> atbaseplot(ring,@plenvelope);
   |  >> atplot(ring,@plenvelope);     (obsolete)
   | 
   |   See also atbaseplot

.. py:function:: CurlyHlindata

   | function [H,Hv]=CurlyHlindata(lindata)
   | 
   |  computes Curly H (dispersion invariant)
   | 
   |  lindata is the first ouptut of [lindata,~,~]=atlinopt(...)
   |  (include at least 2 ouptut arguments for dispersion computation)
   | 
   |  output
   | 
   |  H  : horizontal dispersion invariant
   |  Hv : vertical dispersion invariant
   | 
   | 

.. py:function:: plotbetadisp

   | function [s,plotdata]=plotbetadisp(ring,dpp,plotfun,varargin)
   | PLOTBETADISP Plot beta functions and dispersion
   | 
   | Helper function for tplot:
   | - beta functions on the left axis
   | - diespersion on the right axis
   | 
   |   USAGE:
   |  >> atbaseplot(ring,@plotbetadisp);
   |  >> atplot(ring,@plotbetadisp);      (obsolete)
   | 
   | See also atplot atbaseplot

.. py:function:: plPolynomBComp

   | PLPOLYNOMBCOMP PlotBn coefficient with normalization
   |  PLPOLYNOMBCOMP default plotting function for ATPLOT
   | 
   |  Plots polynomB for ring and ring1

.. py:function:: plClosedOrbit

   | PLCLOSEDORBIT Plots H and V 4x4 closed orbit
   | 
   | Helper function for atplot: plot
   | - H and V closed orbits on left axis
   | - Dispersion on right axis
   | 
   |   EXAMPLEs
   |  >> atbaseplot(ring,@plClosedOrbit,{'dp',0.01});
   |  >> atplot(ring,@plClosedOrbit,'dp',0.01);     (obsolete)
   | 
   |   See also atplot atbaseplot

.. py:function:: plotB0curlyh

   | PLOTB0CURLYH  Plot B and H
   | 
   |  USAGE:
   |  >> atbaseplot(ring,@PLOTB0CURLYH);
   |  >> atplot(ring,@PLOTB0CURLYH);     (obsolete)
   | 
   |   See also atbaseplot

.. py:function:: CurlyH

   | function [H,Hv]=CurlyH(RING,dp,ind)
   | 
   |  computes Curly H (dispersion invariant)
   | 
   |  RING :at lattice
   |  dp  : energy deviation
   |  ind : reference positions
   | 
   |  output
   | 
   |  H  : horizontal dispersion invariant
   |  Hv : vertical dispersion invariant
   | 
   | 

.. py:function:: plBeamSize

   | PLBEAMSIZE Plot H and V beam size
   | 
   |  USAGE:
   |  >> atbaseplot(ring,@PLBEAMSIZE);
   |  >> atplot(ring,@PLBEAMSIZE);     (obsolete)
   | 
   |   See also atbaseplot

.. py:function:: plotWdispP

   | plotWdispP    Plot W functions
   | 
   | Helper function for atplot: plot
   | - W functions (derivatives of beta-functions versus momentum) on left axis
   | - derivative of dispersion on right axis

