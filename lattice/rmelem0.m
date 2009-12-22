function [NEWLATTICE, SHIFTEDINDEX] = rmelem0(LATTICE,ELEMINDEX) 
%RMELEM0 removes elements of length 0 from the accelerator lattice
% NEWLATTICE = RMELEM0(LATTICE,ELEMINDEX)
% [NEWLATTICE, SHIFTEDINDEX] = RMELEM0(LATTICE,ELEMINDEX)
%
% The number of elements in the modified lattice is 
% reduced by length(ELEMINDEX)
%
% SHIFTEDINDEX points to elements in the NEWLATTICE that
% immediateley followed the removed elements in the original LATTICE
%
% See also: SPLITDRIFT MERGEDRIFT

% make a new (empty) lattice
OLDN = length(LATTICE);
NEWN = OLDN-length(ELEMINDEX);

THERING = cell(1,NEWN);
NEWINDEX = zeros(size(ELEMINDEX));

ALLELEMS = 1:length(LATTICE);
ALLELEMS(ELEMINDEX)=0;
NZ = find(ALLELEMS);

for i = 1:length(NZ)
    NEWLATTICE{i}=LATTICE{NZ(i)};
end


% current warning state 
originalwarnstate = warning;
warning off backtrace

    newlength=findspos(NEWLATTICE,NEWN+1);
    oldlength=findspos(LATTICE,OLDN+1);
    if newlength~=oldlength
        error('Elements to be removed must have 0 physical length');
    end

    SHIFTEDINDEX = 1+(ELEMINDEX(:)-(1:length(ELEMINDEX))');
    if size(SHIFTEDINDEX)~=size(ELEMINDEX)
        SHIFTEDINDEX=SHIFTEDINDEX';
    end
    if ELEMINDEX(end)==OLDN
        warning('Last element in the original LATTICE was removed');
    end

    i = 1;
    while i<=length(ELEMINDEX)
        if ~strcmp(LATTICE{ELEMINDEX(i)}.PassMethod,'IdentityPass')
            warning('One or more removed elements used PassMethod other than ''IdentityPass''!');
            i = length(ELEMINDEX)+1;
        else
            i=i+1;
        end
    end
    warning(originalwarnstate);
