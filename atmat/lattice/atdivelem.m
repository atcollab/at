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
if isfield(elem,'T1')
    line(2:end)=cellfun(@(elem) rmfield(elem,'T1'),line(2:end),'UniformOutput',false);
end
if isfield(elem,'R1')
    line(2:end)=cellfun(@(elem) rmfield(elem,'R1'),line(2:end),'UniformOutput',false);
end
if isfield(elem,'T2')
    line(1:end-1)=cellfun(@(elem) rmfield(elem,'T2'),line(1:end-1),'UniformOutput',false);
end
if isfield(elem,'R2')
    line(1:end-1)=cellfun(@(elem) rmfield(elem,'R2'),line(1:end-1),'UniformOutput',false);
end
if isfield(elem,'BendingAngle')
    line=atsetfieldvalues(line,'BendingAngle',elem.BendingAngle*frac(:)/sum(frac));
end
if isfield(elem,'EntranceAngle')
    drentrangle=zeros(size(line));
    drentrangle(1)=elem.EntranceAngle;
    line=atsetfieldvalues(line,'EntranceAngle',drentrangle);
end
if isfield(elem,'ExitAngle')
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
    line=atsetfieldvalues(line,'FringeInt2',fringe2);
end
if isfield(elem,'FringeBendEntrance')
    fringe1=zeros(size(line));
    fringe1(1)=elem.FringeBendEntrance;
    line=atsetfieldvalues(line,'FringeBendEntrance',fringe1);
end
if isfield(elem,'FringeBendExit')
    fringe2=zeros(size(line));
    fringe2(end)=elem.FringeBendExit;
    line=atsetfieldvalues(line,'FringeBendExit',fringe2);
end
if isfield(elem,'FringeQuadEntrance')
    fringe1=zeros(size(line));
    fringe1(1)=elem.FringeQuadEntrance;
    line=atsetfieldvalues(line,'FringeQuadEntrance',fringe1);
end
if isfield(elem,'FringeQuadExit')
    fringe2=zeros(size(line));
    fringe2(end)=elem.FringeQuadExit;
    line=atsetfieldvalues(line,'FringeQuadExit',fringe2);
end

if isfield(elem,'FullGap')
    fringe2=zeros(size(line));
    fringe2([1,end])=elem.FullGap;
    line=atsetfieldvalues(line,'FullGap',fringe2);
end
if isfield(elem,'EdgeEffect1')
    fringe2=zeros(size(line));
    fringe2(1)=elem.EdgeEffect1;
    line=atsetfieldvalues(line,'EdgeEffect1',fringe2);
end
if isfield(elem,'EdgeEffect2')
    fringe2=zeros(size(line));
    fringe2(end)=elem.EdgeEffect2;
    line=atsetfieldvalues(line,'EdgeEffect2',fringe2);
end

end
