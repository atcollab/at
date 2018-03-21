function ring=atApplyVariation(ring,Variables,NewValues)
%ATAPPLYVARIATION Applies a series of variations to variables as described in
% Variables.
%
% Variables is a cell array of structures
%
%
% Variables  struct('Indx',{[indx],...
%                           @(ring,varval)fun(ring,varval,...),...
%                          },...
%                   'Parameter',{{'paramname',{M,N,...},...},...
%                                [initialvarval],...
%                               },...
%                   'LowLim',{[val],[val],...},...
%                   'HighLim',{[val],[val],...},...
%                   )
%
%
% varargin{1} may be a set of delta_0 to apply
%
%
% delata_0: is the applied set of perturbations
% varNum  : length(delta_0)
%

% history of changes
% created 25-8-2012
% updated 28-8-2012 Variables may have a vector as PERTURBINDX
% updated 30-8-2012 added 'Fam' filed for Families
% updated 10-9-2012 Variation.('FUN') Variation.('FUN_PAR') added 
%                   to apply a modification of parameters non esplicitly
%                   stated in the RING structure
%                   Variation.('FUN')=@funtoeval 
%                   must have at least 2 args (THERING,parameter2vary)
% updated 17-9-2012 Variation.('FUN_PAR') function may have more
%                   input parameters that do not change
% updated 17-9-2012 Variation.('FUN_VAR') function may have more
%                   input parameters that may vary
% updated 18-9-2012 Variables(var_indx).('StartVALUE') vector of intial
%                   function variables values
% updated 21-3-2013 main Variables structure change removed familiy option
% updated 23-3-2013 removed loops
% updated 25-3-2013 varibles as absolute values and not variations.
%                   Indx and Parmaeter switched in case of function.
%                   setfield(...Parameter{:}) instead of Parameter{1} or{2}


ok=arrayfun(@setval,Variables,NewValues); %#ok<NASGU>

    function ok=setval(v,val)
        if isa(v.Indx,'function_handle')
            ring=v.Indx(ring,val{1});
        else
            ring(v.Indx)=cellfun(@(elem) setfield(elem,v.Parameter{:},...
                val{1}),ring(v.Indx),'UniformOutput',false);
        end
        ok=0;
    end
        
end

