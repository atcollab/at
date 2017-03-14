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
if length(frac)>1
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
        % Correct bending angles
        line=atsetfieldvalues(line,'BendingAngle',elem.BendingAngle*frac(:)/sum(frac));
        % Keep entrance angle only in first slice
        if isfield(elem,'EntranceAngle')
            drentrangle=zeros(size(line));
            drentrangle(1)=elem.EntranceAngle;
            line=atsetfieldvalues(line,'EntranceAngle',drentrangle);
        end
        % Keep exit angle only on last slice
        if isfield(elem,'ExitAngle')
            drexitangle=zeros(size(line));
            drexitangle(end)=elem.ExitAngle;
            line=atsetfieldvalues(line,'ExitAngle',drexitangle);
        end
        % Suppress dipole fringe effects at intermediate points
        line(2:end)=atsetfieldvalues(line(2:end),'FringeBendEntrance',0);
        line(1:end-1)=atsetfieldvalues(line(1:end-1),'FringeBendExit',0);
    end
    % Suppress quadrupole fringe effects at intermediate points
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
end

end
