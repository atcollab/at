.. _atmatch_module:

atmatch
=======

.. py:module:: atmatch

   ATMATCH

.. rubric:: Functions


.. list-table::

   * - :func:`atApplyVariation`
     - Applies a series of variations to variables as described in
   * - :func:`atDisplayConstraintsChange`
     - This funciton evaluates the contraints defined in Constraints for lattice
   * - :func:`atDisplayVariableChange`
     - this functions retrives variable Values for two rings to compare
   * - :func:`atEvaluateConstraints`
     - This funciton evaluates the contraints defined in Constraints for lattice
   * - :func:`atGetPenalty`
     - 
   * - :func:`atGetVariableValue`
     - 
   * - :func:`atVariableBuilder`
     - atVarBuilder   create a simple variable structure for use with atmatch
   * - :func:`atlinconstraint`
     - Generate a constraint on linear optics for atmatch
   * - :func:`atmatch`
     - function [...

.. py:function:: atApplyVariation

   | Applies a series of variations to variables as described in
   |  Variables.
   
   |  Variables is a cell array of structures
   
   
   |  Variables  struct('Indx',{[indx],...
   |                            @(ring,varval)fun(ring,varval,...),...
   |                           },...
   |                    'Parameter',{{'paramname',{M,N,...},...},...
   |                                 [initialvarval],...
   |                                },...
   |                    'LowLim',{[val],[val],...},...
   |                    'HighLim',{[val],[val],...},...
   |                    )
   
   
   |  varargin{1} may be a set of delta_0 to apply
   
   
   |  delata_0: is the applied set of perturbations
   |  varNum  : length(delta_0)
   

.. py:function:: atDisplayConstraintsChange

   | This funciton evaluates the contraints defined in Constraints for lattice
   |  Constraints: cell array of struct('Fun',@functname,'Min',min,'Max',max,'OtherParameters',otherargs}
   
   |  Constraints: structure array struct('Fun',@(ring)functname(ring,parameters),
   |                                    'Min',min, % of unweigthed funct val
   |                                    'Max',max,
   |                                    'Weight',w,
   |                                    'RefPoints',[])
   
   |  Constraints: structure array struct(...
   |                      'Fun',@(ring,lindata,globaldata,refpts)functname(...
   |                              ring,lindata,globaldata,refpts,parameters),
   |                                    'Min',min, % of unweigthed funct val
   |                                    'Max',max,
   |                                    'Weight',w,
   |                                    'RefPoints',refpts);
   
   |  lindata is the output of atlinopt
   |  globdata.tune=tune fromk atlinopt
   |  globdata.chrom=chrom from atlinopt
   
   |  functname: handle to vector valued function: [res]=functname(THERING,otherargs)
   
   |  min and max have to be the same size as res. (this allows to give functions as limits!)
   
   |  created 30-8-2012
   |  updated 12-10-2012 other_function_args is a cell array, ifit is not it is
   |                     transformed in a cell array

.. py:function:: atDisplayVariableChange

   | this functions retrives variable Values for two rings to compare
   
   |  Variables is a structure array
   
   |  Variables  struct('Indx',{[indx],...
   |                            @(ring,varval)fun(ring,varval,...),...
   |                           },...
   |                    'Parameter',{{'paramname',{M,N,...},...},...
   |                                 [initialvarval],...
   |                                },...
   |                    'LowLim',{[val],[val],...},...
   |                    'HighLim',{[val],[val],...},...
   |                    )
   

.. py:function:: atEvaluateConstraints

   | This funciton evaluates the contraints defined in Constraints for lattice
   |  THERING
   
   |  Constraints: structure array struct(...
   |                      'Fun',@functname(ring,lindata,globaldata),
   |                      'Min',min,
   |                      'Max',max,
   |                      'Weight',w,
   |                      'RefPoints',refpts);
   
   |  lindata is the output of atlinopt at the requested locations
   |  globdata.fractune=tune fromk atlinopt
   |  globdata.chromaticity=chrom from atlinopt
   
   |  functname must return a row vector of values to be optimized
   
   |  min, max and weight must have the same size as the return value of
   |  functname
   

.. py:function:: atGetPenalty

   
   |  Evaluate the penalty function (distance from the target value of every constraint)
   

.. py:function:: atGetVariableValue

   
   |  this functions retrives variable Values
   
   |  Variables is a structure array of structures
   
   
   |  Variables  struct('Indx',{[indx],...
   |                            @(ring,varval)fun(ring,varval,...),...
   |                           },...
   |                    'Parameter',{{'paramname',{M,N,...},...},...
   |                                 [initialvarval],...
   |                                },...
   |                    'LowLim',{[val],[val],...},...
   |                    'HighLim',{[val],[val],...},...
   |                    )
   

.. py:function:: atVariableBuilder(refpts,parameter,highlim,lowlim)

   | atVarBuilder   create a simple variable structure for use with atmatch
   
   |  Single variable : it corresponds to a scalar numeric value to be varied in
   |  the optimization process. It may be applied to several elements.It is
   |  represented as a scalar structure.
   
   |    **var=atVariableBuilder(refpts,parameter,highlim,lowlim)**
   |        refpts:     indices of the variable elements or logical mask
   |        parameter:	cell array defining the field name and indices of the
   |                    variable parameter
   |        lowlim:     minimum parameter value (default: no limit)
   |        highlim:    maximum parameter value (default: no limit)
   
   |        Example:	qf=atgetcells(ring,'FamName','QF');
   |                    **var=atVariableBuilder(qf,{'polynomb',{2}})**;
   
   |    **var=atVariableBuilder(@func,inival,highlim,lowlim)**
   |        func:       function building a new ring for the given variable value
   |                    called as new_ring=func(base_ring,variable)
   |        inival:     initial value of the variable
   |        lowlim:     minimum parameter value (default: no limit)
   |        highlim:    maximum parameter value (default: no limit)
   
   |        Example: **var=atVariableBuilder(@(r,v) some_function(r,v,...), 0.0)**;
   
   |    **var=atVariableBuilder(ring,location,...)**
   |        In this syntax, the location may be specified as the family name of the
   |        variable elements
   
   |        Example: **var=atVariableBuilder(ring,'qf',{'polynomb',{2}})**;
   
   |  Multiple variables: if location,parameter,lowlim and highlim are cell arrays
   |  with the same length or with length 1, **atVariableBuilder** will build a
   |  structure array of variables. Examples:
   
   |    **vars=atVariableBuilder(ring,{'qd','sf'},{{'polynomb',{1,2}},{'polynomb',{1,3}}})**;
   
   |    qf=atgetcells(ring,'FamName','QF');
   |    qd=atgetcells(ring,'FamName','QD');
   |    **vars=atVariableBuilder({qf,qd},{{'polynomb',{1,2}}})**;
   
   |    **vars=atVariableBuilder({qf,@buildring},{{'polynomb',{1,2}},0.0})**
   
   |  More sophisticated variables, can be defined using directly the variable
   |  structure. The general variable definition is:
   
   |  ex: Variab=struct('Indx',{findcells(RING,'FamName','QFM'),...
   |                             k1start(1)},...
   |                    'LowLim',{[],[]},...
   |                    'HighLim',{[],[]},...
   |                    'Parameter',{{'PolynomB',{1,2}},...
   |                                 {'FUN',...
   |                            @(RING,K1Val)VaryQuadFam(RING,K1Val,'QDM')}}...
   |                   );
   

.. py:function:: atlinconstraint(refpts,params,vmin,vmax,weight)

   | Generate a constraint on linear optics for atmatch
   
   | **constraint=atlinconstraint(refpts,params,vmin,vmax,weight)**
   
   | REFPTS Row vector of selected positions
   | PARAMS Cell array describing the desired value at each position
   |        The length of params must be 1 or length(REFPTS)
   |        Each element of PARAMS is itself a cell array defining the field
   |        name and indices in the structure returned by atlinopt. Additional
   |        field names are: 'tune' and 'chromaticity'.
   | VMIN   Minimum value for the constraint
   | VMAX   Maximum value for the constraint
   
   | CONSTRAINT Row structure array to be used in atmatch
   
   |  REFPTS, PARAMS, VMIN, VMAX, WEIGHT must have the same length,
   |        or have length 1
   
   |  Example:
   | >> **c1=atlinconstraint(1,{{'closedorbit',{3}},{'closedorbit',{4}}},[0 0],[0 0],[1/6 6])**;
   
   | See also :func:`atmatch`, :func:`atVariableBuilder`

.. py:function:: atmatch

   | function [...
   |     NewRing,...
   |     penalty,...
   |     dmin...
   |     ]=**atmatch**(...
   |     Ring,...
   |     Variables,...
   |     Constraints,...
   |     Tolerance,...
   |     Calls,...
   |     verbose,...
   |     minimizer,...
   |     twissin)
   
   |  this functions modifies the Variables (parameters in THERING) to obtain
   |  a new THERING with the Constraints verified
   
   |  Ring        : at lattice structure
   |  Variables   : a structure array of parameters to vary with step size.
   |  Constraints : a structure array
   |  Tolerance   : square sum of distance to wished constraints at which the minimizer stops
   |  Calls       : number of calls
   |  verbose     : verbosity 0-3 (see later)
   |  minimizer   : @fminsearch (default) or @lsqnonlin
   |  twissin     : open line matching initial parameters
   |  dpp         : ...,'dpp',0.0,... use atlinopt off energy by dpp
   |  UseParallel : ...,'UseParallel',false,... use parallel pool for
   |                optimization
   
   |  Variables  struct('Indx',{[indx],...
   |                            @(ring,varval)fun(ring,varval,...),...
   |                           },...
   |                    'Parameter',{{'paramname',{M,N,...},...},...
   |                                 [initialvarval],...
   |                                },...
   |                    'LowLim',{[val],[val],...},...
   |                    'HighLim',{[val],[val],...},...
   |                    )
   
   
   |  Constraints: structure array struct(...
   |                      'Fun',@functname(ring,lindata,globaldata),
   |                      'Min',min,
   |                      'Max',max,
   |                      'Weight',w,
   |                      'RefPoints',refpts);
   
   |  lindata is the output of atlinopt at the requested locations
   |  globaldata.fractune=tune from atlinopt
   |  globaldata.chromaticity=chrom from atlinopt
   
   |  functname must return a row vector of values to be optimized
   
   |  min, max and weight must have the same size as the return value of
   |  functname
   
   |  verbose to print out results.
   |                                0 (no output)
   |                                1 (initial values)
   |                                2 (iterations)
   |                                3 (result)
   
   |  Variables are changed within the range min<res<max with a Tolerance Given
   |  by Tolerance
   
   | See also :func:`atlinconstraint`, :func:`atVariableBuilder`

