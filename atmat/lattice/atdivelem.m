function line = atdivelem(elem,frac)
%LINE=ATDIVELEM(ELEM,FRAC) divide an element into pieces
%
%ELEM:  Element to be divided
%FRAC:  Length of each piece, as a fraction of the total length
%
%LINE:  Cell array 
%
% The sum of frac elements does not need to be 1, however for bending
% magnets, the length will be divided according to FRAC, but the total
% bending angle will be divided according to FRAC/SUM(FRAC) so that the
% total bending angle is kept.
%
% Example:
%
%>> qf=atquadrupole('QF',0.1,0.5);
%>> line=atdivelem(qf,[0.5;0.5]); % Split a quadrupole in two halves
%
% See also ATINSERTELEMS ATSLICE ATSPLITELEM

line=atsetfieldvalues(repmat({elem},length(frac),1),'Length',elem.Length*frac(:));
if isfield(elem,'BendingAngle')
    line=atsetfieldvalues(line,'BendingAngle',elem.BendingAngle*frac(:)/sum(frac));
end
if isfield(elem,'EntranceAngle')
    drentrangle=zeros(size(line));
    drentrangle(1)=elem.EntranceAngle;
    line=atsetfieldvalues(line,'EntranceAngle',drentrangle);
end
if isfield(elem,'BendingAngle')
    drexitangle=zeros(size(line));
    drexitangle(end)=elem.ExitAngle;
    line=atsetfieldvalues(line,'ExitAngle',drexitangle);
end

end
