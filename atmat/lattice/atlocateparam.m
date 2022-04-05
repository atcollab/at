function idx = atlocateparam(ring,varargin)
%ATLOCATEPARAM    Private function. Locate the RingParam element
%
% IDX=ATLOCATEPARAM(ring)
%	Extract the RingParam element of the lattice, if present, or create
%	a new one from the lattice elements
%
% ring:         Ring structure
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
    location=find(atgetcells(ring(:,1),'Class','RingParam'), 1);
end
idx = location;
