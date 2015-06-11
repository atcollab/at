function line=atsplitelem(baseelem,varargin)
%ATSPLITELEM Creates a line by inserting one or more elements into a base element
%
%LINE=ATSPLITELEM(BASEELEM,FRAC1,ELEM1[,FRAC2,ELEM2...])
%   Each inserted element is associated with a location given by 0<=FRAC<=1
%   LINE is a cell array containing the sequence of resulting elements
%
% FRACi may be a scalar or a line vector of locations. ELEMi must be a
% single element, or a cell array of elements with the same length as FRACi.
%
% if FRAC = 0, the element is inserted before BASEELEM (no splitting)
% if FRAC = 1, the element is inserted after BASEELEM (no splitting)
%
% if ELEMi = [], nothing is inserted, only the splitting takes place
%
% ATSPLITELEM will split BASEELEM in elements with negative lengths if
% necessary. Elements with length 0 are not generated. When splitting
% dipoles, bending angles are distributed as lengths, and face angles are
% set to zero at splitting points.
%
%
% Examples:
%
%>> dr=atdrift('DR',2);       % Insert quadrupoles inside a drift
%>> qf=atquadrupole('QF',0.1,0.5);
%>> qd=atquadrupole('QD',0.1,-0.5);
%>> line=atsplitelem(dr,0.2,qf,0.8,qd);
%
%>> mk=atmarker('MK');
%>> line=atsplitelem(qf,0,mk);   % Insert a marker before a quadrupole
%
%>> line=atsplitelem(qf,0.5,[]); % Split a quadrupole in two halves
%
% See also ATINSERTELEMS ATDIVELEM

elfrac=varargin(1:2:end);
elems=cell(size(elfrac));
for i=1:length(elfrac)
    ellist=varargin{2*i};
    if ~iscell(ellist)
        ellist={ellist};
    end
    elems{i}(1:length(elfrac{i}))=ellist(:);
end
elems=cat(2,elems{:})';
elfrac=cat(2,elfrac{:})';
ellg=0.5*atgetfieldvalues(elems,'Length')./baseelem.Length;
ellg(isnan(ellg))=0;
drfrac=[elfrac-ellg;1]-[0;elfrac+ellg];
long=drfrac~=0;

drifts=cell(size(drfrac));
drifts(long)=atdivelem(baseelem,drfrac(long));

list=cell(length(drifts)+length(elems),1);
list(1:2:end)=drifts;
list(2:2:end)=elems;

keep=true(size(list));
keep(1:2:end)=long;                         % remove useless elements
keep(2:2:end)=atgetcells(elems,'FamName');  % remove dummy inserts

line=list(keep);
end
