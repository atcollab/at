function elem=atM66(fname,varargin)
%ATM66(FAMNAME,M66,PASSMETHOD)
%   atM66 creates an element that applies an arbitrary matrix m66
%
%FAMNAME	family name
%M66        transfer matrix, defaults to Identity(6)]
%PASSMETHOD	tracking function, defaults to 'Matrix66Pass'
%
%ATM66(FAMNAME,ATSTRUCT,PASSMETHOD)
%   atM66 will generate the matrix by calling FINDM66(ATSTRUCT)
%
%ATSTRUCT   AT structure

[rsrc,m66,method]=decodeatargs({eye(6),'Matrix66Pass'},varargin);
[rsrc,L]=getatarg(rsrc,0.0,'Length');
if isstruct(m66)
    m66=findm66(m66);
end
elem=atbaselem(fname,method,'M66',m66,'Class','Matrix66','Length',L,rsrc{:});
