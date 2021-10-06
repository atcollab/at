function ring = atSetRingProperties(ring,varargin)
%atSetRingProperties	Add or modify properties of the lattice
%
% newring=atSetRingProperties(ring [,key,value])
%   Add or modify the attributes of the RingParam element of the lattice,
%   Insert a new RingParam element if necessary
%
% Standard properties:
%   FamName         Name of the lattice
%   Energy          Ring energy in eV
%   Periodicity     Number of periods to make 2pi phase advance
%   Particle        particle (Particle object)
%
% Optional properties:
%   cavpts          Location of the main cavities (Used by atsetcavity
%

%See also atGetRingProperties

% Fast access if RingParam is the first element, as usual
if isfield(ring{1},'Class') && strcmp(ring{1}.Class, 'RingParam')
    idx=1;
% Otherwise, look around
else
    idx=find(atgetcells(ring(:,1),'Class','RingParam'), 1);
end
if isempty(idx)     % No RingParam element: insert a new one
    s=warning;                          % Save the warning state
    warning('Off','AT:NoRingParam');    % Disable warning
    props=atGetRingProperties(ring);
    warning(s);                         % Restore the warning state
    parmelem=atringparam(props.FamName,props.Energy,props.Periodicity,varargin{:});
    ring=[{parmelem};ring];
else                % Found RingParam, update it
    parms=struct(varargin{:});
    for fn=fieldnames(parms)
        fnn=fn{1};
        ring{idx}.(fnn)=parms.(fnn);
    end
end
end

