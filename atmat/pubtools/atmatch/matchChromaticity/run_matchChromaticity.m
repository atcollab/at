%  match chrmaticity


load dba.mat

variabsxt=atVariableBuilder(RING,{'SF','SD'},{{'PolynomB',{1,3}}});

ConstrChrom=[...
    atlinconstraint(1,{{'chromaticity',{1}}},0,0,1)...
    atlinconstraint(1,{{'chromaticity',{2}}},0,0,1)];

tol=1e-8;
RINGchrom0=atmatch(RING,variabsxt,ConstrChrom,tol,1000,4);

