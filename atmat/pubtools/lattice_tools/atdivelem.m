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
line(2:end)=cellfun(@(elem) rmfield(elem,{'T1','R1'}),line(2:end),'UniformOutput',false);
line(1:end-1)=cellfun(@(elem) rmfield(elem,{'T2','R2'}),line(1:end-1),'UniformOutput',false);
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
if isfield(elem,'FringeInt1')
    fringe1=zeros(size(line));
    fringe1(1)=elem.FringeInt1;
    line=atsetfieldvalues(line,'FringeInt1',fringe1);
end
if isfield(elem,'FringeInt2')
    fringe2=zeros(size(line));
    fringe2(end)=elem.FringeInt2;
    line=atsetfieldvalues(line,'FringeInt1',fringe2);
end

end
