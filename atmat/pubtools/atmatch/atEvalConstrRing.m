function Val=atEvalConstrRing(R,v,c,d)

[R,~,ld0]=atApplyVariation(R,v,d);

[~,Val]=atGetPenalty(R,c);

if ld0>length(Val)
   error('More Variables then Constraints')
end

return