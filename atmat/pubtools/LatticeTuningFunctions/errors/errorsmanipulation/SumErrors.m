function rsum=SumErrors(r1,r2,magindex,indBPM)
% rsum=SumErrors(r1,r2,magindex)
% 
% gets errors from r1 and r2 and sum them into rsum
% 

if length(r1)~=length(r2)
    error('r1 and r2 must be the same lattice with different errors!')
end

if nargin<3
    magindex=1:length(r1);
end

[X1,Y1,S1,T1,R1,P1,bpm1]=GetExistingErrors(r1,magindex);
[X2,Y2,S2,T2,R2,P2,bpm2]=GetExistingErrors(r2,magindex);

bpms.offsetx=bpm1.offsetx+bpm2.offsetx;
bpms.offsety=bpm1.offsety+bpm2.offsety;
bpms.rotation=bpm1.rotation+bpm2.rotation;

rsum=SetExistingError(r1,magindex,indBPM,...
    X1+X2,...
    Y1+Y2,...
    S1+S2,...
    T1+T2,...
    R1+R2,...
    P1+P2,...
    bpms...
    );


return