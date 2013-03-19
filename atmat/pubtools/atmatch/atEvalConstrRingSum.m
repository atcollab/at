function Val=EvalConstrRingSum(R,v,c,d)

R=ApplyVariation(R,v,d);

[~,p]=GetPenalty(R,c);

Val=sum(p);
return