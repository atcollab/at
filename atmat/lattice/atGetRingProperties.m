function varargout = atGetRingProperties(ring,varargin)
%ATGETRINGPROPERTIES Get the ring properties
%
% [V1,V2,...]=ATGETRINGPROPERTIES(RING,'param1','param2',...)
%   Extract lattice properties. Extract from the RingParam element of the
%   lattice if present, or from the lattice elements. The number of outputs
%   corresponds to the number of property names in input.
%
% RING:             Ring structure
%
% Standard properties:
%   'FamName'               Name of the lattice
%   'name'                   "   "   "     "
%   'Energy'                Ring energy [eV]
%   'energy'                 "     "    
%   'Periodicity'           Number of periods to build a full ring
%   'periodicity'             "    "     "    "   "    "  "    "
%   'Particle'              particle (Particle object)
%   'particle'                  "         "      "
%   'cavpts'                Location of the main cavities
%   'beta'                  Relativistic beta of the particles
%   'gamma'                 Relativistic gamma of the particles
%   'rf_frequency'          RF frequency (main cavities) [Hz]
%   'rf_timelag'            RF timelag (main cavities) [m]
%   'BRho'                  Particle rigidity [T.m]
%   'mcf'                   Momentum compaction factor "alpha"
%   'slip_factor'           Slip factor "eta"
%
% Properties for the full ring ('periods' x cells):
%   'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
%   'harmonic_number'          "       "
%   'Circumference'         Ring ircumference [m] (cell_length * Periodicity)
%   'circumference'         "
%   'rf_voltage'            RF voltage [V] (cell_rf_voltage * Periodicity)
%   'revolution_frequency'  Revolution frequency [Hz] (cell_revolution_frequency / Periodicity)
%
% Properties for one cell:
%   'cell_harmnumber'       Harmonic number (cell)
%   'cell_length'           Cell length [m]
%   'cell_rf_voltage'       RF voltage [V] (main cavities)
%   'cell_revolution_frequency' Revolution frequency [Hz] (cell)
%
% Custom properties may be added with atSetRingProperties
%
% Example:
% >> [energy, beta] = atGetRingProperties(ring,'Energy','beta');
%
% PROPERTIES=ATGETRINGPROPERTIES(RING)
%
% RING:             Ring structure
% PROPERTIES:       Structure with fields:
%   'FamName'               Name of the lattice
%   'Energy'                Ring energy [eV]
%   'Periodicity'           Number of periods to build a full ring
%   'Particle'              particle (Particle object)
%   'HarmNumber'            Harmonic number (cell_harmnumber * Periodicity)
%   'cell_harmnumber'       Harmonic number (cell)
%   'cavpts'                Location of the main cavities
%
% For fast access, the ring properties are stored in a RingParam element
% ideally located in the 1st position of the lattice. Without such element,
% the properties are deduced from the lattice contents. This is much slower
% and ATGETRINGPROPERTIES displays a warning indicating how to add the
% RingParam element:
%
%>>ring=atSetRingProperties(ring)
%
%  See also ATSETRINGPROPERTIES

idx=atlocateparam(ring);
if isempty(idx)
%   t1='Slow access to properties because there is no RingParam element.';
%   t2='Consider adding it with the command: ">> ring=atSetRingProperties(ring)".';
%   warning('AT:NoRingParam', '%s\n%s', t1, t2);
    props=struct();
else
    props=rmfield(ring{idx},{'Length','Class','PassMethod'});
end
if isempty(varargin)
    [props,~]=atparamscan(ring,props,'FamName','Energy','Periodicity',...
        'Particle','cell_harmnumber','cavpts');
    props.Particle=atparticle.loadobj(props.Particle);
    props.HarmNumber=props.Periodicity*props.cell_harmnumber;
    varargout={props,idx};
else
    [~,varargout]=atparamscan(ring,props,varargin{:});
end
end
