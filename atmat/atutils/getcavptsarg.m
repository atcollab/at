function varargout = getcavptsarg(args,ring,props)
%GETCAVPTSARG Compute CAVPTS if not provided
%
%[CAVPTS,VARARGS]=GETCAVPTSARG(VARARGS,RING)
%
%VARARGS=GETCAVPTSARG(VARARGS,RING)
%   {'cavpts', CAVPTS} is added to VARARGS
%
% If CAVPTS is not provided, CAVPTS is taken from the RING properties if
% existing, otherwise it uses all cavities

[cavpts,args]=getoption(args,'cavpts',[]);
if isempty(cavpts)
    if nargin < 3, props=atGetRingProperties(ring); end
    try         % Look for cavities in the lattice properties
        cavpts=props.cavpts;
    catch       % Take all cavities
        cavpts=atgetcells(ring,'Frequency');
    end
end

if nargout == 1
    args=[args {'cavpts', cavpts}];
    varargout={args};
else
    varargout={cavpts,args};
end
end