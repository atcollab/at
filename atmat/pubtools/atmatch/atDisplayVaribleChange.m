function atDisplayVaribleChange(THERING1,THERING2,Variables)
% this functions retrives variable Values for two rings to compare
% 
% Variables is a cell array of structures
% 
% Variables  struct('PERTURBINDX',indx,...
%                   'PVALUE',stepsize,
%                   'Fam',boolean,
%                   'FIELD',fieldtochange, 
%                   'IndxInField',{M,N,...})
% 
% created 30-8-2012


% vary varaible
for var_indx=1:length(Variables)
  if ~strcmp(Variables{var_indx}.('FIELD'),'macro') % is not a funcition
      % every variable may have more perturbed indexes.
     PERTURB=Variables{var_indx}.PERTURBINDX;
    
   for i=1:length(PERTURB)
    varname=getcellstruct(THERING1,'FamName',PERTURB(i));
    fieldname=Variables{var_indx}.('FIELD');
    
    value1 = ...
        getfield(THERING1{PERTURB(i)},...
        Variables{var_indx}.('FIELD'),...
        Variables{var_indx}.('IndxInField')...
        );
    
    value2 = ...
        getfield(THERING2{PERTURB(i)},...
        Variables{var_indx}.('FIELD'),...
        Variables{var_indx}.('IndxInField')...
        );
    
    disp([varname{:} '     ' fieldname '     ' num2str(value1) '    ' num2str(value2) '    ' num2str(value2-value1)]);
    
   
   end % end loop variables indexes
  else
      for i=1:length(Variables{var_indx}.('StartVALUE'))  % parameters may be vectors
           
          
         disp([func2str(Variables{var_indx}.('FUN')) '                   ' num2str(i) '     ' num2str(Variables{var_indx}.('StartVALUE')(i)) '    ' num2str(0) '    ' num2str(0)]);
     
      end
  end

end % end loop variables



return
