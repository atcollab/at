function r=VaryMagLength(r,DL,indP,indM)

ndp=length(indP);
ndm=length(indM);

r=setcellstruct(r,'Length',indP,getcellstruct(r,'Length',indP)+DL/ndp);
r=setcellstruct(r,'Length',indM,getcellstruct(r,'Length',indM)-DL/ndm);

return
