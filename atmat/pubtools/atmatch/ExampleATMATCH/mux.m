function [mx]=mux(Seq,indx)
% get value of phase advance for  Seq(indx)


if length(indx)>1
    T=twissring(Seq,0,indx);
    m=cat(1,T.mu);
    mx=m(:,1)/(2*pi);
else
    T=twissring(Seq,0,1:indx+1);
    m=cat(1,T.mu);
    mx=m(end,1)/(2*pi);
end

end
