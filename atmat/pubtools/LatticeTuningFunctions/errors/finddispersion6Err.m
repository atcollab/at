function dispersion=finddispersion6Err(RING, indbpm,indrfc,alpha,delta,inCOD)
%FINDDISPERSION6ERR Gets 6D dispersion with bpm reading errors
%
%  See also findorbit6Err

f0=RING{indrfc(1)}.Frequency;

% plus delta
RINGp=atsetfieldvalues(RING,indrfc,'Frequency',f0-alpha*(+delta)*f0);
orbitp = findorbit6Err(RINGp, indbpm, inCOD);

RINGm=atsetfieldvalues(RING,indrfc,'Frequency',f0-alpha*(-delta)*f0);
orbitm = findorbit6Err(RINGm, indbpm, inCOD);

dispersion=(orbitp-orbitm)/2/delta;


end