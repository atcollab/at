function DiffMat = quantumDiff(elems, varargin)
%QUANTUMDIFF    Compute the radiation-diffusion matrix
%
%DIFFMAT=QUANTUMDIFF(RING)
%   RING:       Closed ring AT structure, containing radiative elements and
%               RF cavity. Radiative elements are identified by a
%               PassMethod ending with 'RadPass'.
%
%DIFFMAT=QUANTUMDIFF(LINE,RADINDEX,ORBITIN)    (Deprecated syntax)
%DIFFMAT=QUANTUMDIFF(...,'orbit',ORBITIN)
%   RADINDEX:   Ignored
%   ORBITIN:    Initial 6-D closed orbit.
%               In this mode, LINE may be a section of a ring.

[orb0,varargs]=getoption(varargin, 'orbit',[]);
[radindex,orb0,varargs]=getargs(varargs,[],orb0); %#ok<ASGLU>

if isempty(orb0)
    orb0=findorbit6(elems);
end

BCUM=atdiffmat(elems,'orbit',orb0);

DiffMat=(BCUM+BCUM')/2;
end
