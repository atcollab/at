function [parmelem, idx] = atfindparam(ring,varargin)
%ATFINDPARAM    Private function. Extract the RingParam element or create a default one
%
% properties=ATFINDPARAM(ring)
%	Extract the RingParam element of the lattice, if present, or create
%	a new one from the lattice elements
%
% ring:         Ring structure
%
%  See also ATGETRINGPROPERTIES, ATSETRINGPROPERTIES

TWO_PI_ERROR = 1.e-4;

idx = atlocateparam(ring);

try
    newparms = struct(varargin{:});
catch
    error('AT:WrongOptions', 'Unexpected options')
end
if isfield(newparms, 'Particle')
    particle = newparms.Particle;
    if ischar(particle) || isstring(particle)
        particle = atparticle(particle);
    end
    newparms.Particle=saveobj(particle);
end

if ~isempty(idx)            % Found RingParam: use it
    parmelem = ring{idx};
    old_nper = parmelem.Periodicity;
    parmelem = strupdate(parmelem, newparms);
    new_nper = parmelem.Periodicity;
    if ~isfield(parmelem, 'Particle')
        parmelem.Particle = saveobj(atparticle('relativistic'));
    end
    if isfield(newparms, 'HarmNumber')
        check_h(newparms.HarmNumber, parmelem.Periodicity);
    elseif ~isfield(parmelem, 'HarmNumber')
        if isfield(parmelem, 'cavpts')
            maincavs = ring(props.cavpts);
        else
            maincavs = findmaincav(ring(atgetcells(ring,'Frequency')));
        end
        if ~isempty(maincavs) && maincavs(1).Frequency ~= 0
            gamma = parmelem.Energy / parmelem.Particle.rest_energy;
            h = round(parmelem.Periodicity * maincavs(1).Frequency / cellfrev(ring, gamma));
            parmelem.HarmNumber = h;
        end
    elseif new_nper ~= old_nper
        parmelem.HarmNumber = parmelem.HarmNumber/old_nper*new_nper;
    end
else                        % No RingParam element : create a new one
%     t1='Slow access to properties because there is no RingParam element.';
%     t2='Consider adding it with the command: ">> ring=atSetRingProperties(ring)".';
%     warning('AT:NoRingParam', '%s\n%s', t1, t2);
    
    cavities = atgetcells(ring(:,1),'Frequency');

    % Look for name
    if isfield(newparms, 'FamName')
        name = newparms.FamName;
        newparms = rmfield(newparms, 'FamName');
    else
        name = '';
    end

    % Look for particle
    if isfield(newparms, 'Particle')
        particle = newparms.Particle;
        newparms = rmfield(newparms, 'Particle');
    else
        particle = saveobj(atparticle('relativistic'));
    end

    % Look for energy
    if isfield(newparms, 'Energy')
        energy = newparms.Energy;
        newparms = rmfield(newparms, 'Energy');
    else
        E0s = atgetfieldvalues(ring,'Energy');
        E0s = E0s(isfinite(E0s));         % Discard undefined values
        maincav=findmaincav(ring(cavities,1));
        if ~isempty(maincav) && isfield(maincav, 'Energy')
            energy=maincav.Energy;
        elseif length(unique(E0s)) == 1
            energy = unique(E0s);
        elseif length(unique(E0s)) > 1
            error('AT:NoEnergy','Energy field not equal for all elements')
        else
            energy=NaN;
        end
    end

    % Look for periodicity
    if isfield(newparms, 'Periodicity')
        nbper = newparms.Periodicity;
        newparms = rmfield(newparms, 'Periodicity');
    else
        % Guess number of periods based on the total bending angle
        dipoles  = atgetcells(ring(:,1),'BendingAngle');
        theta    = atgetfieldvalues(ring(dipoles),'BendingAngle');
        nbp=2*pi/sum(theta);
        nbper=round(nbp);
        if ~isfinite(nbp)
            warning('AT:WrongNumberOfCells','No bending in the cell, ncells set to 1');
            nbper=1;
        elseif abs(nbp-nbper) > TWO_PI_ERROR
            warning('AT:WrongNumberOfCells','non integer number of cells: ncells = %g -> %g',nbp,nbper);
        end
    end

    % Look for harmonic number
    if isfield(newparms, 'HarmNumber')
        h = check_h(newparms.HarmNumber, nbper);
        newparms = rmfield(newparms, 'HarmNumber');
        hargs = {'HarmNumber', h};
    else
        maincav=findmaincav(ring(cavities,1));
        hargs = {};
        if ~isempty(maincav) && maincav(1).Frequency ~= 0
            gamma = energy / particle.rest_energy;
            h = round(nbper * maincav.Frequency/cellfrev(ring, gamma));
            hargs = {'HarmNumber', h};
        end
    end

    % Create the RingParam element
    parmelem=atringparam(name,energy,nbper,'Particle',particle,hargs{:});
    parmelem=strupdate(parmelem,newparms);
end

    function maincav=findmaincav(cavities)
        % Find the main cavity (lowest frequency)
        if isempty(cavities)
            maincav=[];
        else
            [~, ii]=sort(atgetfieldvalues(cavities, 'Frequency'));
            maincav=cavities{ii(1)};
        end
    end

    function frev=cellfrev(ring, gamma)
        % Extract the harmonic number from the main cavity
        vel = sqrt(1-1/gamma/gamma)*PhysConstant.speed_of_light_in_vacuum.value;
        frev = vel / findspos(ring, length(ring)+1);
    end

    function str = strupdate(str, str2)
        % Update a struct with the contents of another one
        f = fieldnames(str2);
        for in=1:length(f)
            fn=f{in};
            str.(fn) = str2.(fn);
        end
    end

    function h = check_h(h, nbper)
        % Check that the harmonic number is a multiple of the periodicity
        if mod(h, nbper) ~= 0
%             error('AT:HarmNumber', 'The harmonic number must be a multiple of %d', nbper);
        end
    end
end
