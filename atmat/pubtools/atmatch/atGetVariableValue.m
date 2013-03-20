function CurrentVal=atGetVariableValue(THERING,Variables,verbose)
% 
% this functions retrives variable Values
% 
% Variables is a cell array of structures
% 
% Variables  struct('PERTURBINDX',indx,...
%                   'PVALUE',stepsize,
%                   'FIELD',fieldtochange, 
%                   'IndxInField', {M,N,...} 
%                    )
% 
% created 25-8-2012
% updated 03-3-2012 to take in account of macro variables! 

CurrentVal=[];

% vary varaible
for var_indx=1:length(Variables)
   if ~strcmp(Variables{var_indx}.FIELD,'macro')
       % every variable may have more perturbed indexes.
    PERTURB=Variables{var_indx}.PERTURBINDX;
    
    
   value=zeros(1,length(PERTURB));
   
   for i=1:length(PERTURB)
    varname=getcellstruct(THERING,'FamName',PERTURB(i));
    value(i) = ...
            getfield(THERING{PERTURB(i)},...
            Variables{var_indx}.FIELD,...
            Variables{var_indx}.IndxInField...
            );
        
   if verbose
       disp([ varname '  '  Variables{var_indx}.FIELD ' : ' num2str(value(i))]);
   end  
   
   end % end loop variables indexes
   CurrentVal{var_indx}.Val=value;
   
   else
   CurrentVal{var_indx}.Val=Variables{var_indx}.StartVALUE;
   
   end
end % end loop variables



return
