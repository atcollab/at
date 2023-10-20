function idx = atlocateparam(ring)
%ATLOCATEPARAM    Private function. Locate the RingParam element
%
% IDX=ATLOCATEPARAM(ring)
%
% ring:         Ring structure
% IDX:          Index of the 1st RingParam element in the ring
%
%  See also ATGETRINGPROPERTIES, ATSETRINGPROPERTIES

persistent location         % Location saved for fast access

% Assume RingParam in 1st position in the ring
if isempty(location)
    location = 1;
end

% Check if it is where expected, otherwise look for it
if ~(length(ring) >= location && ...
     isfield(ring{location},'Class') && ...
     strcmp(ring{location}.Class, 'RingParam'))
    location=find(atgetcells(ring,'Class','RingParam'), 1);
end
idx = location;
