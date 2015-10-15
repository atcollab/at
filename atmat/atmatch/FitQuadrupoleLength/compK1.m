function dif=compK1(arc,indQ1,indQ2,varargin)

KQ1=getcellstruct(arc,'PolynomB',indQ1,1,2);
KQ2=getcellstruct(arc,'PolynomB',indQ2,1,2);

dif=KQ2(1)-KQ1(1);

return