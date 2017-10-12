function line = atdivelem(elem,frac,varargin)
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
% Optional arguments:
% 'KeepAxis', if present, rotations translations are kept at all slices 
%
% See also ATINSERTELEMS ATSLICE ATSPLITELEM ENTRANCEFIELDS EXITFIELDS

% get fields names to keep at entrance and exit only
entfield=entrancefields(varargin{:});
exfield=exitfields(varargin{:});

[entrancef,el]=mvfield(struct(),elem,entfield); % Extract entrance fields
[exitf,el]=mvfield(struct(),el,exfield);           % extract exit fields
                                        % split the bare element
line=atsetfieldvalues(repmat({el},length(frac),1),'Length',el.Length*frac(:));
if isfield(elem,'BendingAngle')
    line=atsetfieldvalues(line,'BendingAngle',el.BendingAngle*frac(:)/sum(frac));
    line=atsetfieldvalues(line,'EntranceAngle',0.0);
    line=atsetfieldvalues(line,'ExitAngle',0.0);
end

line{1}=mvfield(line{1},entrancef);     % Set back entrance fields
line{end}=mvfield(line{end},exitf);     % Set back exit fields
end
