function [parmelem,pvals] = atparamscan(ring,parmelem,varargin)
%ATPARAMSCAN    Private function. Updates the RingParam element
%
%NEWPARMS=ATPARAMSCAN(RING,PARMS,VARARGIN)
%
%  See also ATGETRINGPROPERTIES, ATSETRINGPROPERTIES

global GLOBVAL %#ok<GVMIS>

TWO_PI_ERROR = 1.e-4;

store=struct();

pvals=cellfun(@getany, varargin, 'UniformOutput', false);

    function v=getany(param)
        switch lower(param)
            case {'famname','name'}
                v=getname(ring);
            case 'particle'
                v=getparticle(ring);
            case 'energy'
                v=getenergy(ring);
            case 'periodicity'
                v=getnp(ring);
            case 'cell_harmnumber'
                v=getcell_h(ring);
            case 'cavpts'
                v=getmaincav(ring);
            case 'beta'
                v=get_beta(ring);
            case 'gamma'
                v=get_gamma(ring);
            case 'cell_length'
                v=get_cell_length(ring);
            case {'harmnumber','harmonic_number'}
                v=getnp(ring)*getcell_h(ring);
            case 'circumference'
                v=getnp(ring)*get_cell_length(ring);
            case 'rf_frequency'
                v=get_unique(ring(getmaincav(ring)), 'Frequency');
            case 'cell_rf_voltage'
                v=get_sum(ring(getmaincav(ring)), 'Voltage');
            case 'rf_voltage'
                v=getnp(ring) * get_sum(ring(getmaincav(ring)), 'Voltage');
            case 'rf_timelag'
                v=get_unique(ring(getmaincav(ring)), 'TimeLag');
            case 'cell_revolution_frequency'
                v=get_cell_frev(ring);
            case 'revolution_frequency'
                v=get_cell_frev(ring) / getnp(ring);
            case 'brho'
                clight=PhysConstant.speed_of_light_in_vacuum.value;
                energy=getenergy(ring);
                rest_energy=getparticle(ring).rest_energy;
                if rest_energy == 0
                    rest_energy = 1.0e6*PhysConstant.electron_mass_energy_equivalent_in_MeV.value;
                end
                v=sqrt(energy^2 - rest_energy^2)/clight;
            case 'mcf'
                v=get_mcf(ring);
            case 'slip_factor'
                gamma=get_gamma(ring);
                v=1/gamma/gamma - get_mcf(ring);
            case 'radiation'
                v=getradcav(ring,{'RadPass', 'QuantDiffPass', 'CavityPass'});
            case 'active_cavity'
                v=getradcav(ring,'CavityPass');
            otherwise
                v=parmelem.(param);
        end
    end

    function found=getradcav(ring, pattern)
        % disp('compute radon, cavon');
        found=false;
        for i=1:length(ring)
            passmethod=ring{i}.PassMethod;
            if endsWith(passmethod, pattern)
                found=true;
                break;
            end
        end
    end

    function alphac=get_mcf(ring)
        if isfield(store,'mcf')
            alphac=store.mcf;
        else
            [~,rnorad]=check_radiation(ring,false,'force');
            alphac=mcf(rnorad);
            store.mcf=alphac;
        end
    end

    function l=get_cell_length(ring)
        if isfield(store,'cell_length')
            l=store.cell_length;
        else
            % disp('Compute length');
            l=findspos(ring,length(ring)+1);
            store.cell_length=l;
        end
    end

    function frev=get_cell_frev(ring)
        if isfield(store,'cell_revolution_frequency')
            frev=store.cell_revolution_frequency;
        else
            % disp('Compute frev');
            vel = get_beta(ring)*PhysConstant.speed_of_light_in_vacuum.value;
            frev = vel / get_cell_length(ring);
            store.cell_revolution_frequency=frev;
        end
    end

    function name=getname(~)
        % Look for name
        if isfield(parmelem,'FamName')
            name=parmelem.FamName;
        else
            % disp('Compute name');
            if isfield(GLOBVAL,'LatticeFile')
                name = GLOBVAL.LatticeFile;
            else
                name = '';
            end
            parmelem.FamName=name;
        end
    end

    function energy=getenergy(ring)
        % Look for energy
        if isfield(parmelem,'Energy')
            energy=parmelem.Energy;
        else
            % disp('Compute energy');
            cavities=ring(getmaincav(ring));
            if ~isempty(cavities) && isfield(cavities{1}, 'Energy')
                energy=cavities{1}.Energy;
            elseif isfield(GLOBVAL,'E0')
                energy=GLOBVAL.E0;
            else
                E0s = atgetfieldvalues(ring,'Energy');
                E0s = E0s(isfinite(E0s));   % Discard undefined values
                E0s = unique(E0s);
                if length(E0s) == 1
                    energy = E0s;
                elseif length(E0s) > 1
                    error('AT:NoEnergy','Energy field not equal for all elements')
                else
                    energy=NaN;
                end
            end
            parmelem.Energy=energy;
        end
    end        

    function cell_h=getcell_h(ring)
        % Look for ring harmonic number
        if isfield(parmelem,'HarmNumber')
            [parmelem, ring_h]=strpop(parmelem,'HarmNumber');
        end
        % Look for cell harmonic number
        if isfield(parmelem,'cell_harmnumber')
            cell_h=parmelem.cell_harmnumber;
        elseif exist('ring_h','var')
            cell_h = ring_h/getnp(ring);
            parmelem.cell_harmnumber=cell_h;
        else
            maincav=ring(getmaincav(ring));
            if ~isempty(maincav)
                if maincav{1}.Frequency ~= 0
                    cell_h = round(maincav{1}.Frequency/get_cell_frev(ring));
                elseif isfield(maincav{1},'HarmNumber')
                    cell_h = maincav{1}.HarmNumber/getnp(ring);
                else
                    cell_h=NaN;
                end
            else
                cell_h = NaN;
            end
            parmelem.cell_harmnumber=cell_h;
        end
    end

    function maincav=getmaincav(ring)
        if isfield(parmelem,'cavpts')
            maincav=parmelem.cavpts;
        else
            % disp('Compute maincav');
            maincav=atmaincavities(ring);
            parmelem.cavpts=maincav;
        end
    end

    function nbper=getnp(ring)
        % Look for periodicity
        if isfield(parmelem,'Periodicity')
            nbper=parmelem.Periodicity;
        else
            % Guess number of periods based on the total bending angle
            % disp('Compute periodicity');
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
            parmelem.Periodicity=nbper;
        end
    end

    function particle=getparticle(~)
        % Look for particle
        if isfield(parmelem,'Particle')
            partstruct = parmelem.Particle;
        else
            % disp('Compute particle');
            partstruct = saveobj(atparticle('relativistic'));
            parmelem.Particle=partstruct;
        end
        particle=atparticle.loadobj(partstruct);
    end

    function gamma=get_gamma(ring)
        if isfield(store,'gamma')
            gamma=store.gamma;
        else
            rest_energy=getparticle(ring).rest_energy;
            if rest_energy == 0
                rest_energy = 1.0e6*PhysConstant.electron_mass_energy_equivalent_in_MeV.value;
            end
            gamma=getenergy(ring) / rest_energy;
            store.gamma=gamma;
        end
    end

    function beta=get_beta(ring)
        rest_energy=getparticle(ring).rest_energy;
        gammainv=rest_energy / getenergy(ring);
        beta= sqrt(1.0 - gammainv/gammainv);
    end

    function [str,varargout]=strpop(str,varargin)
        varargout=cellfun(@extr,varargin,'UniformOutput',false);
        function v=extr(nm)
            v=str.(nm);
            str=rmfield(str,nm);
        end
    end

    function value=get_unique(cavities,attrname)
        if isempty(cavities)
            value=NaN;
        else
            value=unique(atgetfieldvalues(cavities,attrname));
            if length(value) > 1
                error('AT:Unique','The attribute "%s" is different among the selected cavities',attrname);
            end
        end
    end

    function value=get_sum(cavities,attrname)
        if isempty(cavities)
            value=NaN;
        else
            value=sum(atgetfieldvalues(cavities,attrname));
        end
    end
end
