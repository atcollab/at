function elem=atEnergyRamp(fname,varargin)
%ATENERGYRAMP Creates an energy ramp element with Class 'EnergyRamp'
%ATENERGYRAMP(FAMNAME,TurnsRamp,EnergyRamp)
%
%  INPUTS
%  1. FAMNAME	   - Family name
%  2. TurnsRamp    - array (1,NPointsRamp) with turn numbers with change of
%                    slope for theenergy ramp
%  3. EnergyRamp   - Energy at the turn TurnsRamp
%
%  OUTPUTS
%  1. ELEM - Structure with the AT element
%
%  EXAMPLES
%  1. atEnergyRamp('EnRamp',[0,100,1000],[200e6,220e6,1000e6]);
%


[rsrc,turnsramp, energyramp] = decodeatargs({0,200e6},varargin);
[method,rsrc]     = getoption(rsrc,'PassMethod','EnergyRampPass');
[cl,rsrc]         = getoption(rsrc,'Class','EnergyRamp');

NPointsRamp=min([length(turnsramp),length(energyramp)]);

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',0,...
    'TurnsRamp',turnsramp,'EnergyRamp',energyramp,...
    'NPointsRamp',NPointsRamp,rsrc{:});

end
