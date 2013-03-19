function [ConstVal,penalty,check]=atGetPenalty(R,Constraints)
% 
% Evaluate Constraints values and returns their value, and the
% penalty function (distance from the target value of every constraint)
% and the check vector, a boolean vetor of converged (1) and not converged
% (0) constraints
% returns also upper and lower bound. 
% The bounds should be of the same lenght as the output of 'Fun'
%
%
% created 25(?)-8-2012
% updated 28-12-2012 (updated EvaluateConstraints to handle vectors.)

[ConstVal,LowLim,UpLim]=atEvaluateConstraints(R,Constraints);


check=zeros(size(ConstVal));
penalty=zeros(size(ConstVal));

% disp([LowLim' ConstVal' UpLim'])

% check constraints and get penalty
for i=1:length(ConstVal)
    if ConstVal(i)<UpLim(i) && ConstVal(i)>LowLim(i)
        check(i)=1;
    else
        %        disp([LowLim(i) ConstVal(i) UpLim(i)])
        
        if ConstVal(i)>UpLim(i)
            penalty(i)=sqrt((ConstVal(i)-UpLim(i))^2);
        else
            penalty(i)=0;
        end
        if ConstVal(i)<LowLim(i)
            penalty(i)=penalty(i)+sqrt((ConstVal(i)-LowLim(i))^2);
        else
            penalty(i)=0+penalty(i);
        end
        
        check(i)=0;
    end
end



