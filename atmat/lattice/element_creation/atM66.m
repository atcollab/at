function elem=atM66(fname,varargin)
%ATM66  Create an element applying an arbitrary 6x6 transfer matrix
%
%ATM66(FAMNAME, M66, M66RAD, PASSMETHOD)
%
%FAMNAME	family name
%M66        transfer matrix, defaults to Identity(6)
%M66RAD     transfer matrix when radiation is activated. Defaults to M66
%PASSMETHOD	tracking function, defaults to 'Matrix66Pass'
%
%ATM66(FAMNAME,ATSTRUCT,PASSMETHOD)
%   atM66 will generate the matrix by calling FINDM66(ATSTRUCT)
%
%ATSTRUCT   AT structure

[rsrc,m66,m66rad,method] = decodeatargs({eye(6),[],'Matrix66Pass'},varargin);
[method,rsrc]     = getoption(rsrc,'PassMethod',method);
[cl,rsrc]         = getoption(rsrc,'Class','Matrix66');

if isstruct(m66)
    m66=findm66(m66);
end
if isempty(m66rad)
    m66rad=m66;
end

elem=atbaselem(fname,method,'Class',cl,'M66',m66,'M66Rad',m66rad,rsrc{:});

end
