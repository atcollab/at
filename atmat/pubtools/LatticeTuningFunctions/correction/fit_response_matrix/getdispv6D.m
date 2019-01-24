function dv=getdispv6D(r,ib,indrfc,alpha,delta,inCOD)

d=finddispersion6Err(r,ib,indrfc,alpha,delta,inCOD);

f0=r{indrfc(1)}.Frequency;

dv=-d(3,:)'./(alpha*f0) ; % [m/Hz]

end