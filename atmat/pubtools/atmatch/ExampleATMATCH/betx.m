function [bx]=betx(Seq,indx)
% get value of betx for  Seq(indx)

T=twissring(Seq,0,indx);
b=cat(1,T.beta);
bx=b(:,1)';

end