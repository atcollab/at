function [THERING,delta_0,varNum,lowlim,highlim]=atApplyVariation(THERING,Variables,varargin)
% this functions applys a series of varaiation to variables as described in
% Variables.
%
% Variables is a cell array of structures
%
% Variables  struct('PERTURBINDEX',indx,...
%                   'PVALUE',stepsize,       % starting Perturbation value
%                   'StartVALUE',start, 
%                   'Fam',1 or 0, % family
%                   'FIELD',fieldtochange,
%                   'FUN',@funtoeval,        % return THERING given by this function.
%                   'FUN_PAR',{parameters},  % return THERING given by this function.
%                   'IndInField',{M,N,...}
%                    )
%
% varargin{1} may be a set of delta_0 to apply
%
%
% delata_0: is the applied set of perturbations
% varNum  : length(delta_0)
%
% created 25-8-2012
% updated 28-8-2012 Variables may have a vector as PERTURBINDX
% updated 30-8-2010 added 'Fam' filed for Families
% updated 10-9-2010 Variation.('FUN') Variation.('FUN_PAR') added 
%                   to apply a modification of parameters non esplicitly
%                   stated in the RING structure
%                   Variation.('FUN')=@funtoeval 
%                   must have at least 2 args (THERING,parameter2vary)
% updated 17-9-2010 Variation.('FUN_PAR') function may have more
%                   input parameters that do not change
% updated 17-9-2010 Variation.('FUN_VAR') function may have more
%                   input parameters that may vary
% updated 18-9-2010 Variables{var_indx}.('StartVALUE') vector of intial
%                   function variables values


vi=1;
delta_0=[];
lowlim=[];
highlim=[];
   
% vary varaible
for var_indx=1:length(Variables)
    
     if ~strcmp(Variables{var_indx}.('FIELD'),'macro') %isempty(Variables{var_indx}.('FUN')) % is not a funcition
      % every variable may have more perturbed indexes.
      PERTURB=Variables{var_indx}.('PERTURBINDX');
     
  
        for i=1:length(PERTURB)
            if length(varargin)==1
                variation=varargin{1}(vi);
                
            else
                variation=Variables{var_indx}.('PVALUE');
                % count variables
            end
            
            
            
            delta_0(vi)=variation; %#ok<AGROW>
            if isfield(Variables{var_indx},'LowLim') 
                if ~isempty(Variables{var_indx}.('LowLim'))
                lowlim(vi)=Variables{var_indx}.('LowLim'); %#ok<AGROW>
                end
            else
                lowlim(vi)=[];
            end
            if isfield(Variables{var_indx},'HighLim')
                if ~isempty(Variables{var_indx}.('HighLim'))
                highlim(vi)=Variables{var_indx}.('HighLim'); %#ok<AGROW>
                end
            else
                highlim(vi)=[];
            end
            
            if ~Variables{var_indx}.('Fam')
                vi=vi+1;
            end
            
           oldvalue = ...
                getfield(THERING{PERTURB(i)},...
                Variables{var_indx}.('FIELD'),...
                Variables{var_indx}.('IndxInField')...
                );
            
            
            THERING{PERTURB(i)} = ...
                setfield(THERING{PERTURB(i)},...
                Variables{var_indx}.('FIELD'),...
                Variables{var_indx}.('IndxInField'),...
                oldvalue+variation);
            
            
            
        end % end loop variables indexes
        if Variables{var_indx}.('Fam')
            vi=vi+1;
        end
    else % is a function
        
        % function may take a vector as input variation parameters.
        %nparam=Variables{var_indx}.('FUN_VAR');
        %funinputvariables={}; %zeros(length(nparam),1);
        %length(nparam) %  number of cells.
        %for par=1:1 %length(nparam)
        vecvar=[];
        for i=1:length(Variables{var_indx}.('StartVALUE'))  % parameters may be vectors
            
            if length(varargin)==1
                variation=varargin{1}(vi);
                delta_0(vi)=variation; %#ok<AGROW>
                
            else
                variation=Variables{var_indx}.('PVALUE');
                
                delta_0(vi)=Variables{var_indx}.('StartVALUE')(i)+variation; %#ok<AGROW>
            end
            
            
            if isfield(Variables{var_indx},'LowLim')
                if ~isempty(Variables{var_indx}.('LowLim'))
                    lowlim(vi)=Variables{var_indx}.('LowLim'); %#ok<AGROW>
                end
            else
                lowlim(vi)=[];
            end
            if isfield(Variables{var_indx},'HighLim')
                if ~isempty(Variables{var_indx}.('HighLim'))
                    highlim(vi)=Variables{var_indx}.('HighLim'); %#ok<AGROW>
                end
            else
                highlim(vi)=[];
            end
            
            %vecvar=[vecvar variation]; %#ok<AGROW>
            vecvar=[vecvar, delta_0(vi)]; %#ok<AGROW>
            
            vi=vi+1;
            
        end
        
        %funinputvariables=[funinputvariables vecvar]; %#ok<AGROW>
        
        function_handle=Variables{var_indx}.('FUN');
        
        THERING=feval(function_handle,THERING,vecvar);
        
        if Variables{var_indx}.('Fam')
            vi=vi+1;
        end
    end
end % end loop variables

varNum=length(delta_0);

return
