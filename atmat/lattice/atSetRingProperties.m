function ring = atSetRingProperties(ring,varargin)
%atSetRingProperties	Add or modify properties of the lattice
%
% newring=atSetRingProperties(ring [,key,value]...)
%   Add or modify the attributes of the RingParam element of the lattice,
%   Insert a new RingParam element if necessary
%
% Available properties (property names are case-independent):
%   'FamName'               Name of the lattice
%   'name'                   "   "   "     "
%   'Energy'                Ring energy [eV]
%   'Periodicity'           Number of periods to build a full ring
%   'Particle'              particle (perticle name or Particle object)
%   'cavpts'                Location of the main cavities
%   'rf_frequency'          RF frequency (main cavities) [Hz]. Use 'nominal'
%                           to set the nominal frequency
%   'rf_timelag'            RF timelag (main cavities) [m]
%
% Properties for the full ring ('periods' x cells):
%   'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
%   'harmonic_number'          "       "
%   'rf_voltage'            RF voltage [V] (cell_rf_voltage * Periodicity)
%
% Properties for one cell:
%   'cell_harmnumber'       Harmonic number (for 1 cell)
%   'cell_rf_voltage'       RF voltage per cell [V] (main cavities)
%
% Additional custom fields may be added. They can be retrieved by
% atGetRingProperties and are saved in files. Custom field names are
% case-dependent.
%
% For fast access, some ring properties are stored in a RingParam element
% ideally located in the 1st position of the lattice. If there is no such
% element, atSetRingProperties will add it.
%
% See also atGetRingProperties

parms=struct();
rfparms=struct();
hparms=struct();
cellfun(@canonical,varargin(1:2:end),varargin(2:2:end),'UniformOutput',false);

idx = atlocateparam(ring);
if isempty(idx)
    % No RingParam element: create and insert a new one
    idx=1;
    parmelem=struct('FamName','','Length',0,'PassMethod','IdentityPass','Class','RingParam');
    parmelem=strupdate(parmelem,parms);
    ring=[{parmelem};ring];
    % Add missing properties
    [parmelem,~]=atparamscan(ring,parmelem,...
        'Energy','Periodicity','Particle','cell_harmnumber','cavpts');
else
    % Update the existing RingParam element
    parmelem = ring{idx};
    [parmelem] = strupdate(parmelem, parms);
    % Add possibly missing properties
    [parmelem,~]=atparamscan(ring,parmelem,'Particle','cell_harmnumber','cavpts');
end

% Delay the setting of harmonic number after the periodicity is known
if isfield(hparms, 'HarmNumber')
    nperiods=parmelem.Periodicity;
    parmelem.cell_harmnumber=check_h(nperiods,hparms.HarmNumber);
elseif isfield(hparms, 'cell_harmnumber')
    nperiods=parmelem.Periodicity;
    parmelem.cell_harmnumber=check_h(nperiods,nperiods*hparms.cell_harmnumber);
end
ring{idx}=parmelem;

rfparms=expand(rfparms);
if ~isempty(rfparms)
    ring=atsetcavity(ring,'cavpts',parmelem.cavpts,rfparms{:});
end

    function str = strupdate(str, str2)
        % Update a struct with the contents of another one
        f = fieldnames(str2);
        for in=1:length(f)
            name=f{in};
            str.(name) = str2.(name);
        end
    end

    function cell_h=check_h(nper,ring_h)
        % Check the validity of the harmonic number
        cell_h=ring_h/nper;
        % Check on full ring
        if isfinite(ring_h) && (round(ring_h) - ring_h) ~= 0
            error('AT:Invalid','The harmonic number must be integer')
        end
%       % Check on cell
%       if (round(cell_h) - cell_h) ~= 0
%           error('AT:Invalid','The harmonic number must be a multiple of %i',nper)
%       end
    end

    function canonical(pname,value)
        if ~(ischar(pname) || isstring(pname))
            error('AT:WrongOptions', 'Unexpected option')
        end
        switch lower(pname)
            case {'length','class','beta','gamma','cell_length','circumference',...
                    'cell_revolution_frequency','revolution_frequency','brho',...
                    'mcf','slip_factor','radiation','is_6d','has_cavity'}
                error('AT:Invalid','The property "%s" is read-only',pname);
            case {'famname','name'}
                parms.FamName=value;
            case 'energy'
                parms.Energy=value;
            case 'periodicity'
                parms.Periodicity=check_periodicity(value);
            case 'particle'
                parms.Particle=check_particle(value);
            case {'harmnumber','harmonic_number'}
                hparms.HarmNumber=value;
            case 'cell_harmnumber'
                hparms.cell_harmnumber=value;
            case 'rf_frequency'
                rfparms.Frequency=value;
            case 'rf_voltage'
                rfparms.Voltage=value;
            case 'rf_timelag'
                rfparms.TimeLag=value;
            case 'cell_rf_voltage'
                rfparms.cell_voltage=value;
            otherwise
                parms.(pname)=value;
        end
    end

    function vcell=expand(vstruct)
        vcell=reshape([fieldnames(vstruct)';struct2cell(vstruct)'],1,[]);
    end

    function partstruct=check_particle(particle)
        if ischar(particle) || isstring(particle)
            particle = atparticle(particle);
        elseif ~isa(particle,'atparticle')
            error('AT:TypeError', 'Cannot convert argument to particle');
        end
        % Convert Particle object to struct for saving in .m or .mat files
        partstruct=saveobj(particle);
    end


    function nper=check_periodicity(value)
        if ~(isnumeric(value) && value > 0 && value == floor(value))
            error('AT:ValueError','Periodicity must a positive integer');
        end
        bnd=atgetfieldvalues(ring,'BendingAngle');
        totalangle=value*sum(bnd(isfinite(bnd)));
        if totalangle > 2*pi
            warning('AT:ValueWarning',['The total angle of the lattice (%g)',...
                'is larger then 2 pi'], totalangle);
        end
        nper=value;
    end
end
