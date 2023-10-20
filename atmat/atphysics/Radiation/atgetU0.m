function U0 = atgetU0(ring,varargin)
%ATGETU0 Computes Energy loss per turn in eV .
%
%U0=ATGETU0(RING)   Return the energy loss/turn in eV for the full ring.
%
% RING:     Ring structure
% U0:       Energy loss per turn in eV
%
%U0=ATGETU0(...,'periods',PERIODS) Select the number of periods
%
% PERIODS if the number of periods to take into account (Default: full ring)
%
%U0=ATGETU0(...,'method',METHOD)	Choose the method
%
% METHOD:   'integral': (default) The losses are obtained from
%                       Losses = Cgamma / 2pi * EGeV^4 * I2
%                       Takes into account bending magnets and wigglers.
%           'tracking': The losses are obtained by  tracking without cavities.
%                       Needs radiation ON, takes into account all radiating elements.
%
% See also RINGPARA ATSETCAVITY ATENERGY

[energy, nper]=atGetRingProperties(ring,'Energy','Periodicity');
[nper,varargs]=getoption(varargin,'periods',nper);
[method, varargs]=getoption(varargs,'method','integral'); %#ok<ASGLU>
if strcmpi(method, 'integral')
    U0=nper*integral(ring);
elseif strcmpi(method, 'tracking')
    U0=nper*tracking(ring);
else
    error('AT:InvalidMethod', 'Invalid method: %s', method);
end

    function U0=integral(ring)
        % Losses = Cgamma/2/pi*EGeV^4*I2
        e_mass=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;	% eV
        cspeed = PhysConstant.speed_of_light_in_vacuum.value;                   % m/s
        e_radius=PhysConstant.classical_electron_radius.value;                  % m
        Cgamma=4.0*pi*e_radius/3/e_mass^3;                                      % [m/eV^3]
        Brho=sqrt(energy*energy - e_mass*e_mass)/cspeed;
        coef=Cgamma/2/pi*energy^4;
        
        % Dipole radiation
        dipoles  = atgetcells(ring,'BendingAngle');
        theta    = atgetfieldvalues(ring(dipoles),'BendingAngle');
        lendp    = atgetfieldvalues(ring(dipoles),'Length');
        I2d=sum(abs(theta.*theta./lendp));
        % Wiggler radiation
        iswiggler=@(elem) isfield(elem,'Class') && strcmp(elem.Class,'Wiggler') ...
            && ~strcmp(elem.PassMethod,'DriftPass');
        wigglers=cellfun(iswiggler, ring);
        I2w=sum(cellfun(@wiggler_i2,ring(wigglers)));
        % EnergyLoss radiation
        iseloss=@(elem) isfield(elem,'Class') && strcmp(elem.Class,'EnergyLoss') ...
            && ~strcmp(elem.PassMethod,'IdentityPass');
        eloss=cellfun(iseloss, ring);
        I2e=sum(cellfun(@eloss_i2, ring(eloss)));
        % Additional radiation
        I2x=sum(atgetfieldvalues(ring(atgetcells(ring,'I2')),'I2'));
        I2=I2d+I2w+I2e+I2x;            % [m-1]
        U0=coef*I2;                    % [eV]
        
        function i2=wiggler_i2(elem)
            rhoinv=elem.Bmax/Brho;
            coefh=elem.By(2,:)*rhoinv;
            coefv=elem.Bx(2,:)*rhoinv;
            i2=elem.Length*(coefh*coefh'+coefv*coefv')/2;
        end

        function i2=eloss_i2(elem)
            i2=elem.EnergyLoss/coef;
        end
    end

    function U0=tracking(ring)
        % Ensure 6d is enabled
        check_6d(ring,true);
        % Turn cavities off
        ringtmp=atdisable_6d(ring,'allpass','','cavipass','auto');
        o0=zeros(6,1);
        o6=ringpass(ringtmp,o0);
        U0=-o6(5)*energy;
    end

end
