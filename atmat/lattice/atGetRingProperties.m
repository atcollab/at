function props = atGetRingProperties(ring)
%ATGETRINGPROPERTIES Get the ring properties
%
% properties=ATGETRINGPROPERTIES(ring)
%	Extract data from the RingParam element of the lattice, if present,
%	or from the lattice elements
%
% ring:         Ring structure
%
% properties:   Structure with fields:
%   FamName         Name of the lattice
%   Energy          Ring energy in eV
%   Periodicity     Number of periods to make 2pi phase advance
%   Particle        particle (Particle object)
%
%  See also ATSETRINGPROPERTIES

global GLOBVAL
TWO_PI_ERROR = 1.e-4;

% Fast access if RingParam is the first element, as usual
if isfield(ring{1},'Class') && strcmp(ring{1}.Class, 'RingParam')
    idx=1;
    % Otherwise, look around
else
    idx=find(atgetcells(ring(:,1),'Class','RingParam'), 1);
end
if ~isempty(idx)            % Found RingParam
    parmelem=ring{idx};
%     if ~isfield(parmelem, 'Particle')
%         parmelem.Particle=Particle('electron');
%     end
else                        % No RingParam element
    t1='Slow access to properties because there is no RingParam element.';
    t2='Consider adding it with the command: ">> ring=atSetRingProperties(ring)".';
    warning('AT:NoRingParam', '%s\n%s', t1, t2);

    cavities = atgetcells(ring(:,1),'Frequency');
    E0s = atgetfieldvalues(ring,'Energy');
    E0s = E0s(isfinite(E0s));         % Discard undefined values
    if any(cavities) && isfield(ring{find(cavities,1)},'Energy')
        energy=ring{find(cavities,1)}.Energy;
    elseif isfield(GLOBVAL,'E0')
        energy=GLOBVAL.E0;
    elseif length(unique(E0s)) == 1
        energy = unique(E0s);
    elseif length(unique(E0s)) > 1
        error('AT:NoEnergy','Energy field not equal for all elements')
    else
        error('AT:NoEnergy',...
            ['Energy not defined (searched in '...
            '''RingParam'',''RFCavity'',GLOBVAL.E0,',...
            ' field ''Energy'' of each element)']);
    end
    
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
    parmelem=atringparam('',energy,nbper);
end
props=rmfield(parmelem,{'Length','Class','PassMethod'});
end
