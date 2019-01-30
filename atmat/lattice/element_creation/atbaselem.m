function elem = atbaselem(famname,method,varargin)
%ATBASELEM  Create an AT element structure + various checks
%
%ELEM=ATBASELEM(FAMNAME,METHOD,'FIELDNAME1',VALUE1,...) create AT element
%   Create an AT element structure and check the consistence of
%   PolynomA, PolynomB, MaxOrder and NumIntSteps
%
%  NOTES
%    1. length of PolynomA and PolynomB are equal (zero padding)
%    2. MaxOrder is always lenght(PolynomA) - 1

DefaultNumIntSteps = 10;

[famname,rsrc] = getoption(varargin,'FamName',famname);
[method,rsrc]  = getoption(rsrc,'PassMethod',method);
[lg,rsrc]      = getoption(rsrc,'Length',0);
[defmax,rsrc]  = getoption(rsrc,'DefaultMaxOrder',0);
elem           = struct('FamName',famname,'PassMethod',method,'Length',lg,rsrc{:});

% Making PolynomA of same length with zero padding when necesssary
ab = isfield(elem,{'PolynomA','PolynomB'});
if any(ab)
    if ~ab(1), elem.PolynomA=[]; end
    if ~ab(2), elem.PolynomB=[]; end
    if ~isfield(elem,'MaxOrder')
        elem.MaxOrder=max([defmax+1 find(abs(elem.PolynomB)>0,1,'last') find(abs(elem.PolynomA)>0,1,'last')])-1;
    end
    la = length(elem.PolynomA);
    lb = length(elem.PolynomB);
    if la < elem.MaxOrder+1, elem.PolynomA=[elem.PolynomA zeros(1,elem.MaxOrder+1-la)]; end
    if lb < elem.MaxOrder+1, elem.PolynomB=[elem.PolynomB zeros(1,elem.MaxOrder+1-lb)]; end
    la = length(elem.PolynomA);
    lb = length(elem.PolynomB);
    if la < lb, elem.PolynomA = [elem.PolynomA zeros(1,lb-la)]; end
    if lb < la, elem.PolynomB = [elem.PolynomB zeros(1,la-lb)]; end
    if elem.Length ~= 0 && ~isfield(elem,'NumIntSteps')
        elem.NumIntSteps = DefaultNumIntSteps;
    end
end

end
