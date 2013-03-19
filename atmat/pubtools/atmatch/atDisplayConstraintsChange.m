function atDisplayConstraintsChange(THERING1,THERING2,Constraints)
% This funciton evaluates the contraints defined in Constraints for lattice
% Constraints: cell array of struct('Fun',@functname,'Min',min,'Max',max,'OtherParameters',otherargs}
% 
% 
% functname: handle to vector valued function: [res]=functname(THERING,otherargs)
%
% min and max have to be the same size as res. (this allows to give functions as limits!)
%
% created 30-8-2012
% updated 12-10-2012 other_function_args is a cell array, ifit is not it is
%                    transformed in a cell array


low=[]; % bolean vector
high=[]; % bolean vector

% evaluate constraints and get penalty function.



for constr_indx=1:length(Constraints) % loop contraints to get value and limits.
    ConstrName='                     ';
    
    function_handle=Constraints{constr_indx}.Fun;
    fName=func2str(function_handle);
    
    if length(fName)<length(ConstrName)
        ConstrName(1:length(fName))=fName;
    else
        ConstrName=fName(1:length(ConstrName));
    end
    
    % CstrVal is a vector, output of function_handle.
    CstrVal1=feval(function_handle,THERING1);
    CstrVal2=feval(function_handle,THERING2);

    [~,penalty1]=atGetPenaltyDif(THERING1,Constraints(constr_indx));
    [~,penalty2]=atGetPenaltyDif(THERING2,Constraints(constr_indx));

    low=Constraints{constr_indx}.('Min'); % bolean vector
    high=Constraints{constr_indx}.('Max'); % bolean vector
% %      
%      size(CstrVal1')
%      size(CstrVal2')
%      size(low')
%      size(high')
%      size(penalty1')
%      size(penalty2')
% %     
%     size(repmat(ConstrName,length(CstrVal1),1))
%     size(repmat('  ',length(CstrVal1),1))
%     size(num2str((1:length(CstrVal1))','%d'))
%     size(num2str([ CstrVal1' CstrVal2' low' high' penalty1' penalty2' ],'%03.3e\t'))
%     

% use isrow to determine dimension ununbiguously.
    disp([repmat(ConstrName,length(CstrVal1),1) repmat('  ',length(CstrVal1),1)...
        num2str((1:length(CstrVal1))','%d')  repmat('  ',length(CstrVal1),1) ...
        num2str([ CstrVal1' CstrVal2' low' high' penalty1' penalty2' ],'%03.3e\t')]);
    
    
end
