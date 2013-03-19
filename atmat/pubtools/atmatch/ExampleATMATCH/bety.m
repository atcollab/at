function [by]=bety(Seq,indx)
% get value of bety for  Seq(indx)

T=twissring(Seq,0,indx);
b=cat(1,T.beta);
by=b(:,2)';

end