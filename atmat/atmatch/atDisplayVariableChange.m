function atDisplayVariableChange(ring1,ring2,Variables)
% this functions retrives variable Values for two rings to compare
%
% Variables is a structure array
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
% created 30-8-2012
% updated 21-03-2012 Varaibles are structure array and not cell array of
%                    structures
%                   Indx and Parmaeter switched in case of function.
%                   setfield(...Parameter{:}) instead of Parameter{1} or{2}

disp('Final variable values:')
disp('   ')
disp('Name       field          before    after   variation')
ok=arrayfun(@(v) dispv(v,ring1,ring2),Variables); %#ok<NASGU>
disp('   ')
disp('-----oooooo----oooooo----oooooo----')
disp('   ')

    function ok=dispv(Variable,ring1,ring2)
        if isa(Variable.Indx,'function_handle')
            funcname=sprintf('%-20.20s',func2str(Variable.Indx));
            for i=1:length(Variable.Parameter)
                fprintf('%-23.23s %8g    %8g    %8g\n',...
                    funcname,Variable.Parameter(i),0,0);
            end
        else
            ok=cellfun(@(elem1,elem2) dd(elem1,elem2,Variable.Parameter),...
                ring1(Variable.Indx),ring2(Variable.Indx)); %#ok<NASGU>
        end
        ok=0;
    end

    function ok=dd(elem1,elem2,Parameter)
        value1 = getfield(elem1,Parameter{:});
        value2 = getfield(elem2,Parameter{:});
        fprintf('%-10.10s %-12.12s %8g    %8g    %8g\n',...
            elem1.FamName,Parameter{1},value1,value2,(value2-value1));
        ok=0;
    end
end
