function CurrentVal=atGetVariableValue(ring,Variables)
%
% this functions retrives variable Values
%
% Variables is a structure array of structures
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

% history of changes
% created 25-8-2012
% updated 03-3-2012 to take in account of macro variables!
% updated 25-3-2013 Indx and Parameter switched in case of function.
%                   getfield(...Parameter{:}) instead of Parameter{1} or{2}


CurrentVal=arrayfun(@getval,Variables,'UniformOutput',false);

    function value=getval(v)
        if isa(v.Indx,'function_handle')
            value=v.Parameter;
        else
            value=getfield(ring{v.Indx(1)},v.Parameter{:});
%           value=mean(cellfun(@(elem) getfield(elem,v.Parameter{:}),ring(v.Indx)));
        end
    end

end
