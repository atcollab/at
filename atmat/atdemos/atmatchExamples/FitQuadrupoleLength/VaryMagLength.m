function r=VaryMagLength(r,DL,indP,indM)

ndp=length(indP);
ndm=length(indM);

r=atsetfieldvalues(r,indP,'Length',atgetfieldvalues(r,indP,'Length')+DL/ndp);
r=atsetfieldvalues(r,indM,'Length',atgetfieldvalues(r,indM,'Length')-DL/ndm);

return
