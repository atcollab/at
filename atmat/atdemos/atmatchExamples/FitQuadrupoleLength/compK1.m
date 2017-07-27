function dif=compK1(arc,indQ1,indQ2,varargin)

KQ1=atgetfieldvalues(arc,indQ1,'PolynomB',{1,2});
KQ2=atgetfieldvalues(arc,indQ2,'PolynomB',{1,2});

dif=KQ2(1)-KQ1(1);

return