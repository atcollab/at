function sVal=atEvalConstrRingDifS(R,v,c,d)
% Get sum of squares of distances from target

try
    [R,~,ld0]=atApplyVariation(R,v,d);
    
    %[ConstVal,LowLim,UpLim]=EvaluateConstraints(R,c);
    %Val=min(abs(ConstVal-LowLim),abs(ConstVal-UpLim));

    %     disp('Current Set: ')
    %     disp(d')
    
    [~,Val]=atGetPenaltyDif(R,c);
catch
    disp('This variation is too  large. Applying in 20 steps.')
    try % try to do it in two steps 
        % if a rematch occurs in the Variation this may ease the convergence.
        % it 1
        [~,d0,~]=atApplyVariation(R,v); % get starting values.
        Nst=50;
        for step=1:Nst
        [R,~,ld0]=atApplyVariation(R,v,d0+(d-d0)/Nst*step);
        
        %[ConstVal,LowLim,UpLim]=EvaluateConstraints(R,c);
        %Val=min(abs(ConstVal-LowLim),abs(ConstVal-UpLim));
        disp(['Current Set ' num2str(step) '/' num2str(Nst) ': '])
        disp((d0+(d-d0)/Nst*step)')
        [~,Val]=atGetPenaltyDif(R,c);
        end
        
    catch
        error('Variation is too large, AT may  not determine optics')
    end
end
%
% if ld0>length(Val)
%    error(['More Variables ' num2str(ld0) ' then Constraints ' num2str(length(Val)) ''])
% end
sVal=sum(Val.^2);
%disp(['Constr Val= :' num2str(Val)])
%disp(['sum square= :' num2str(sVal)])
return