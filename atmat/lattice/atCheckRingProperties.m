function [props,idx] = atCheckRingProperties(ring)
%ATCHECKRINGPROPERTIES Get the ring properties if existing
%
% PROPERTIES=ATCHECKRINGPROPERTIES(ring)
%	Extract data from the RingParam element of the lattice, if present,
%	otherwise return an empty structure
%
%  See also ATGETRINGPROPERTIES, ATSETRINGPROPERTIES

idx = atlocateparam(ring);
if ~isempty(idx)
    parmelem=ring{idx};
    if isfield(parmelem, 'Particle')
        particle=atparticle.loadobj(parmelem.Particle);
        props=rmfield(parmelem,{'Length','Class','PassMethod','Particle'});
        props.Particle=particle;
    else
        props=rmfield(parmelem,{'Length','Class','PassMethod'});
    end
else
    props=struct();
end