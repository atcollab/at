function dv=getdispv6D(r,ib,indrfc,alpha,delta,inCOD)

d=finddispersion6Err(r,ib,indrfc,alpha,delta,inCOD);

dv=d(3,:)';

end