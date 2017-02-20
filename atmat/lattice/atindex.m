function indexstruct = atindex(lattice)
%ATINDEX extracts the information about element families and
% indexes from AT lattice
%
%ATI=ATINDEX(LATTICE)
%  returns a srtucture with fields named after element family
%  containing an array of their AT indexes;
%
%   ATI.QF = [...]
%   ATI.QD = [...];
%   ...
%
%ATI=ATINDEX() Uses the global variable THERING
%
% See also ATGETCELLS

if nargin < 1
    indexstruct=ati(THERING);
else
    indexstruct=ati(lattice);
end

    function ids=ati(lattice)
        ids=struct();
        for i=1:length(lattice);
            try
                famname=lattice{i}.FamName;
                try
                    a.(famname)=0;
                catch %#ok<CTCH>
                    try
                        famname=['x' famname]; %#ok<AGROW>
                        a.(famname)=0;
                    catch
                        famname='badname';
                    end
                end
            catch
                famname='unnamed';
            end
            if isfield(ids,famname)
                ids.(famname)(end+1)=i;
            else
                ids.(famname)(1)=i;
            end
        end
    end
end
