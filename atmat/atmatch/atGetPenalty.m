function penalty=atGetPenalty(ConstVal,Constraints)
% 
% Evaluate the penalty function (distance from the target value of every constraint)
%
vals=cat(2,ConstVal{:});
low=cat(2,Constraints.Min);
high=cat(2,Constraints.Max);
weight=cat(2,Constraints.Weight);
penalty=zeros(size(vals));

toohigh=vals>high;
toolow=vals<low;

penalty(toohigh)=(vals(toohigh)-high(toohigh))./weight(toohigh);
penalty(toolow)=(low(toolow)-vals(toolow))./weight(toolow);

% check=~(toohigh|toolow);
end
