function [THERINGout,penalty,dmin]=atmatch(...
    THERING,Variables,Constraints,Tolerance,Calls,algo,verbose)
% this functions modifies the Variables (parameters in THERING) to obtain
% a new THERING with the Constraints verified
%
% Variables   : a cell array of structures of parameters to vary with step size.
% Constraints : a cell array of structures
%
% Variables:   struct('PERTURBINDX',indx,
%                     'PVALUE',stepsize,
%                     'Fam',1,...      % 1=change PERTURBINDX all equal
%                                        0=change PERTURBINDX all different
%                     'LowLim',[],...  % only for lsqnonlin
%                     'HighLim',[],... % only for lsqnonlin
%                     'FIELD',fieldtochange, 
%                     'IndxInField',{{M,N,...}}) 
%
% Variables:   struct('FUN',@(ring,varVAL)functname(ring,varVAL,...),
%                     'PVALUE',stepsize,
%                     'FIELD','macro', 
%                     'StartVALUE',valstart, % size of varVAL  
%                     'IndxInField',{{M,N,...}}) 
% 
% Constraints: struct('Fun',@(ring)functname(ring,...),
%                     'Min',min,
%                     'Max',max,
%                     'Weight',1)
%
% verbose (boolean) to print out end results. 
% 
% Variables are changed within the range min<res<max with a Tolerance Given
% by Tolerance
%
% using least square.
%
%

%
% created : 27-8-2012
% updated : 28-8-2012 constraints 'Fun' may output vectors
% updated : 17-9-2012 output dmin
% updated : 6-11-2012 added simulated annealing (annealing) 
%                     and simplex-simulated-annealing (simpsa)
%                     for global minimum search.
% updated : 21-02-2013 TipicalX included in  
% updated (major) : 11-03-2013 
%                   function named 'atmatch'.
%                   anonimous functions constraints and variable
%                   Variable limits in variable.
%


% get length of variations and higher and lower limits
[~,delta_0,~,Blow,Bhigh]=atApplyVariation(THERING,Variables);


% tolfun is te precisin of the minimum value tolx the accuracy of the
% parameters (delta_0)
% tipicalx is the value of tipical change of a variable 
% (useful if varibles have different ranges)
if isempty(findcells(Variables,'FIELD','macro'))
    tipx=ones(size(delta_0));
    initval=atGetVariableValue(THERING,Variables,0); % get intial variable values
    for ivar=1:length(initval)
        a=initval{ivar};
        if a.Val(1)~=0
            tipx(ivar)=a.Val(1);
        end
    end
    
else
    tipx=ones(size(delta_0)); % the default
end


options = optimset(...
    'Display','iter',...% 
    'MaxFunEvals',Calls*100,...
    'MaxIter',Calls,...
    ...'TypicalX',tipx,...
    'TolFun',Tolerance,...);%,...'Algorithm',{''},...
    'TolX',Tolerance);%,...  %                         
     %'Algorithm',{'levenberg-marquardt',1e-6}

   
switch algo
    case 'lsqnonlin'
        options = optimset('TypicalX',tipx);
% Difference between Target constraints and actual value.
        f = @(d)atEvalConstrRingDif(THERING,Variables,Constraints,d); % vector
    case {'fminsearch','annealing'}
        fs = @(d)atEvalConstrRingDifS(THERING,Variables,Constraints,d); % scalar (sum of squares of f)
end

[~,penalty]=atGetPenaltyDif(THERING,Constraints);

disp(['f²: ' num2str(penalty.^2)]);
disp(['Sum of f²: ' num2str(sum(penalty.^2))]);

%% Least Squares
if sum(penalty)>Tolerance
    % minimize sum(f_i²)
    switch algo
        case 'lsqnonlin'
            
            dmin=lsqnonlin(f,delta_0,Blow,Bhigh,options);
            % dmin=lsqnonlin(f,delta_0,[],[],options);
        
        case 'fminsearch'
            %options = optimset('OutputFcn', @stopFminsearchAtTOL);
            
            dmin = fminsearch(fs,delta_0,options); % wants  a scalar
        
    end
else
    dmin=delta_0;
end
%%

THERINGout=atApplyVariation(THERING,Variables,dmin);
[~,penalty]=atGetPenaltyDif(THERINGout,Constraints);

if verbose
disp('-----oooooo----oooooo----oooooo----')
disp('   ')
disp(['f²: ' num2str(penalty.^2)]);
disp(['Sum of f²: ' num2str(sum(penalty.^2))]);
disp('   ')
disp('-----oooooo----oooooo----oooooo----')
disp('   ')
disp('Final constraints values:')
disp('   ')
disp('Name          lat_indx      before         after           low            high       min dif before    min dif after  ')
atDisplayConstraintsChange(THERING,THERINGout,Constraints);
disp('   ')
disp('-----oooooo----oooooo----oooooo----')
disp('    ')
disp('Final variable values:')
disp('   ')
disp('Name      field      before    after   variation')
atDisplayVaribleChange(THERING,THERINGout,Variables);
disp('   ')
disp('-----oooooo----oooooo----oooooo----')
disp('   ')
end

%atGetVariableValue(THERINGout,Variables,1);
%format long
%disp(dmin')

return