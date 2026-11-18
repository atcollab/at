function rout = elempass(elem,rin,varargin)
%ELEMPASS Tracks particles through a single element
%
%ROUT=ELEMPASS(ELEM,RIN) Tracks particles through ELEM
%
%   ELEM:       lattice element
%   RIN         6xN matrix: input coordinates of N particles
%
%   ROUT        6xN matrix: output coordinates of N particles
%
%  ROUT=ELEMPASS(...,'PassMethod',PASSMETHOD,...)
%     Use PASSMETHOD (default: ELEM.PassMethod)
%
%  ROUT=ELEMPASS(...,'Energy',ENERGY,...)
%     Use ENERGY and ignore the 'Energy' field of elements
%
%  ROUT=ELEMPASS(...,'Particle',PARTICLE,...)
%     Use PARTICLE (default: relativistic)
%
% See also: RINGPASS, LINEPASS

[props.Energy,varargs]=getoption(varargin,'Energy',0.0);
[particle,varargs]=getoption(varargs,'Particle',[]);
[methodname,varargs]=getoption(varargs,'PassMethod',elem.PassMethod); %#ok<ASGLU>
if ~isempty(particle)
    props.Particle=particle;
end

method_req_energy = { ...
    'BeamLoadingCavityPass' , ...
    'CrabCavityPass', ...
    'GWigSymplecticPass', ...
    'RFCavityPass' ...
    };

if endsWith(methodname, {'RadPass', 'QuantPass'}) || any(strcmp(methodname, method_req_energy))
    if props.Energy <= 0
        error("Energy parameter must be defined.");
    end
end

rout = feval(methodname,elem,rin,props);

end