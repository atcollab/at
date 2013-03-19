function [Val,low,high]=atEvaluateConstraints(R,Constraints)
% This funciton evaluates the contraints defined in Constraints for lattice
% THERING
% Penalties : DISTANCE from wished values
% Constraints: cell array of struct('Fun',@(ring)functname(ring,parameters),
%                                   'Min',min, % of unweigthed funct val
%                                   'Max',max,
%                                   'Weight',w,
%                                   )
% 
% 
% functname: handle to vector valued function: [res]=functname(THERING,otherargs)
%
% min and max have to be the same size as res. (this allows to give functions as limits!)
%
% created 26-8-2012
% updated 28-8-2012 vector res handling.
% updated 28-8-2012 add Weigth field.
% updated 11-10-2012 any number of additonal parameters in cell array format.
%                    if only one parameter still ok.
% updated 11-03-2013 proper use of anonytmous functions introduced.


Val=[]; % bolean vector
low=[]; % bolean vector
high=[]; % bolean vector

% check if there is a weigth otherwise set to 1.

% evaluate constraints and get penalty function.
valuesindx=1;
for constr_indx=1:length(Constraints) % loop contraints to get value and limits.
    
    % !!!! CHANGED from MATCH_AT_GN
    
    function_handle=Constraints{constr_indx}.('Fun');
     
    % CstrVal is a vector, output of function_handle.
    CstrVal=feval(function_handle,R); 
       
    % check if there is a weigth otherwise set to 1.
    if not(isfield(Constraints{constr_indx},'Weight'))
        Constraints{constr_indx}.('Weight')=1;
    end
    
    % assign values
    Val(valuesindx:valuesindx+length(CstrVal)-1)=CstrVal.*Constraints{constr_indx}.('Weight');
    
    % assigns values to limits of the constraint
    try
        low(valuesindx:valuesindx+length(CstrVal)-1)=Constraints{constr_indx}.('Min').*Constraints{constr_indx}.('Weight'); % bolean vector
        high(valuesindx:valuesindx+length(CstrVal)-1)=Constraints{constr_indx}.('Max').*Constraints{constr_indx}.('Weight'); % bolean vector
    catch %#ok<CTCH>
        error('size of Constraints{i}.Min and Constraints{i}.Max must be the same of the output of Constraints{i}.Fun')
    end
    
    % update values counter
    valuesindx=valuesindx+length(CstrVal);
    
end
