function sVal=EvalConstrRingDifSimpsa(d,R,v,c)
% Get sum of squares of distances from target
% simpsa function requires d to be the first argument.

try
    [R,~,ld0]=ApplyVariation(R,v,d);
    
    %[ConstVal,LowLim,UpLim]=EvaluateConstraints(R,c);
    %Val=min(abs(ConstVal-LowLim),abs(ConstVal-UpLim));

    %     disp('Current Set: ')
    %     disp(d')
    
    [~,Val]=GetPenaltyDif(R,c);
catch
    disp('This variation is too  large. Applying in 20 steps.')
    try % try to do it in two steps 
        % if a rematch occurs in the Variation this may ease the convergence.
        % it 1
        [~,d0,~]=ApplyVariation(R,v); % get starting values.
        Nst=50;
        for step=1:Nst
        [R,~,ld0]=ApplyVariation(R,v,d0+(d-d0)/Nst*step);
        
        %[ConstVal,LowLim,UpLim]=EvaluateConstraints(R,c);
        %Val=min(abs(ConstVal-LowLim),abs(ConstVal-UpLim));
        disp(['Current Set ' num2str(step) '/' num2str(Nst) ': '])
        disp((d0+(d-d0)/Nst*step)')
        [~,Val]=GetPenaltyDif(R,c);
        end
        
    catch
        error('Variation is too large, AT may  not determine optics')
    end
end
%
% if ld0>length(Val)
%    error(['More Variables ' num2str(ld0) ' then Constraints ' num2str(length(Val)) ''])
% end

sVal=sum(Val.^2)

return