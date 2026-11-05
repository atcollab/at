function element = atswapkickmapdata(element,kickmapKey)
%ATSWAPKICKMAPDATA Activate a stored insertion-device kick map for atinsertiondevicekickmap element
%
% ELEMENT=ATSWAPKICKMAPDATA(ELEMENT,KICKMAPKEY) copies the tracking data from
% ELEMENT.KickmapStore.(KICKMAPKEY) to ELEMENT and sets ActiveKickmap to
% KICKMAPKEY. The PassMethod and all unrelated element fields are unchanged.
% KICKMAPKEY is a string key to stored kickmap data. 
%
% The selected kick map must contain Filename_in, Normalization_energy,
% Nslice, Length, xkick, ykick, xkick1, ykick1, xtable and ytable.
%
% Example:
%   ring{idindex} = atswapkickmapdata(ring{idindex},'LV');
% See also atinsertiondevicekickmap

if ~isstruct(element) || ~isfield(element,'KickmapStore') || ...
        ~isstruct(element.KickmapStore)
    error('AT:atswapkickmapdata:MissingStore', ...
        'The element does not contain a valid KickmapStore.');
end

if ~(ischar(kickmapKey) || (isstring(kickmapKey) && isscalar(kickmapKey)))
    error('AT:atswapkickmapdata:InvalidKey', ...
        'The kick-map key must be a character vector or string scalar.');
end
kickmapKey = char(kickmapKey);

if ~isfield(element.KickmapStore,kickmapKey)
    availableKeys = fieldnames(element.KickmapStore);
    error('AT:atswapkickmapdata:UnknownKey', ...
        'KickmapStore does not contain key "%s". Available keys: %s', ...
        kickmapKey,strjoin(availableKeys,', '));
end

kickmapData = element.KickmapStore.(kickmapKey);
trackingFields = {'Filename_in','Normalization_energy','Nslice','Length', ...
    'xkick','ykick','xkick1','ykick1','xtable','ytable'};
missingFields = trackingFields(~isfield(kickmapData,trackingFields));
if ~isempty(missingFields)
    error('AT:atswapkickmapdata:IncompleteKickmap', ...
        'KickmapStore.%s is missing fields: %s', ...
        kickmapKey,strjoin(missingFields,', '));
end

for field = trackingFields
    element.(field{1}) = kickmapData.(field{1});
end
element.ActiveKickmap = kickmapKey;
end
