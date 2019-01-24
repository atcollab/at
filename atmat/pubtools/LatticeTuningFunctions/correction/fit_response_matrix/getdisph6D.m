function dh=getdisph6D(r,ib,indrfc,alpha,delta,inCOD)

d=finddispersion6Err(r,ib,indrfc,alpha,delta,inCOD);

f0=r{indrfc(1)}.Frequency;

dh=-d(1,:)'./(f0*alpha) ;  % [m/Hz] *Hz

end