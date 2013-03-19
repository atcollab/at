function [ConstVal,penalty,check,LowLim,UpLim]=atGetPenaltyDif(R,Constraints)
% 
% Evaluate Cosntraints values and returns their value, and the
% penalty function (distance from the target value of every constraint)
% and the check vector, a boolean vetor of converged (1) and not converged
% (0) constraints
%

[ConstVal,LowLim,UpLim]=atEvaluateConstraints(R,Constraints);

check=zeros(size(ConstVal));
penalty=zeros(size(ConstVal));

% disp([LowLim' ConstVal' UpLim'])
LowLim(LowLim==0)=-1e-11;
UpLim(UpLim==0)=1e-11;
% check constraints and get penalty
for i=1:length(ConstVal)
    if UpLim(i)==LowLim(i)
        
        penalty(i)=abs(min((ConstVal(i)-LowLim(i)),(ConstVal(i)-UpLim(i))));
     
      % penalty(i)=abs(min((ConstVal(i)-LowLim(i))/LowLim(i),(ConstVal(i)-UpLim(i))/UpLim(i)));
 
    elseif ConstVal(i)<UpLim(i) && ConstVal(i)>LowLim(i)
        penalty(i)=0;
        check(i)=1;
        %disp([LowLim(i) ConstVal(i) UpLim(i) penalty(i)]);
    elseif ConstVal(i)>UpLim(i)
        penalty(i)=abs((ConstVal(i)-LowLim(i))); % too high distance from lower bound
       % penalty(i)=abs((ConstVal(i)-LowLim(i))/LowLim(i)); % too high distance from lower bound
 
        check(i)=0;
    elseif ConstVal(i)<LowLim(i)
        penalty(i)=abs((ConstVal(i)-UpLim(i)));  % too low distance from higher bound
       % penalty(i)=abs((ConstVal(i)-UpLim(i))/UpLim(i));  % too low distance from higher bound
  
        check(i)=0;
    end
    
end


