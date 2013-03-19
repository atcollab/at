function Val=atEvalConstrRingDif(R,v,c,d)


[R,~,ld0]=atApplyVariation(R,v,d);

%[ConstVal,LowLim,UpLim]=EvaluateConstraints(R,c);
%Val=min(abs(ConstVal-LowLim),abs(ConstVal-UpLim));

[~,Val]=atGetPenaltyDif(R,c);

% if ld0>length(Val)
%    error(['More Variables ' num2str(ld0) ' then Constraints ' num2str(length(Val)) ''])
% end

return