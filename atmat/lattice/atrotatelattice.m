function newring = atrotatelattice(ring,index)
%ATROTATELATTICE    Circularly shift the lattice elements
%
%NEWRING = ATROTATELATTICE(RING,INDEX)
%   Return a new lattice such that RING(INDEX) is the first element

if index>length(ring)
    error('index too big!');
end

if index==1
    newring=ring;
else
    ploc = atgetcells(ring,'Class','RingParam');
    newploc=[ploc(index:end);ploc(1:index-1)];
    newring=[ring(index:end);ring(1:index-1)];
    parms = ring(ploc);             % Extract RingParam elements
    if length(parms) > 1            % Keep only the 1st one
        parms = parms(1:1);
    end
    if ~isempty(parms)              % Update the cavpts argument
        if isfield(parms{1}, 'cavpts')
            cavpts = parms{1}.cavpts;
            if islogical(cavpts)
                npts = cavpts;
            else
                npts = false(size(ring));
                npts(cavpts)=true;
            end
            newcavpts=[npts(index:end);npts(1:index-1)];
            parms{1}.cavpts = find([false;newcavpts(~newploc)]);
        end
    end
    % Put back the RingParam element in the 1st position
    newring=[parms;newring(~newploc)];
end