function dh=getdisph6D(r,ib,indrfc,alpha,delta,inCOD)

d=finddispersion6Err(r,ib,indrfc,alpha,delta,inCOD);

dh=d(1,:)';

end