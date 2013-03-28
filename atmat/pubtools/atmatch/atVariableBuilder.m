function variable=atVariableBuilder(r,familyname,fieldtochange,LowLim,HighLim)
% 
% defines a simple variable structure for use with atmatch
% 
% set limits if not given and formats structure correctly.
% 
% familyname={'QD','SF','Dipole1',...}
% fieldtochange={{'PolynomB',{1,2}},{'PolynomB',{1,3}},{'BendingAngle',{1}},...}
% 
% use to create standard variables that change AT lattice structure fields.
% 
% to create more sophisticated variables, please use directly the variable
% structure. The general variable definition is:
%
% ex: Variab=struct('Indx',{findcells(RING,'FamName','QFM'),...
%                            k1start(1)},...
%                   'LowLim',{[],[]},...
%                   'HighLim',{[],[]},...
%                   'Parameter',{{'PolynomB',{1,2}},...
%                                {'FUN',...
%                           @(RING,K1Val)VaryQuadFam(RING,K1Val,'QDM')}}...
%                  );
% 
% if running a matching routine many times, 
% please consider again the definition of variables above. 
% Providing the Indx as an input rether than recalculating it, is faster.
% see example on matching.
% 
% ex : v=atVariableBuilder(r,familyname,fieldtochange)
% ex : v=atVariableBuilder(r,familyname,fieldtochange,LowLim,HighLim)
% ex : variable=atVariableBuilder(RING,...
%                                 {'QD','B'},...
%                                 {{'PolynomB',{1,2}},...
%                                  {'BendingAngle',{1,1}}
%                                 })
%                                 

% history of changes
% created 21-03-2013

variable=[];

for nstru=1:length(familyname)
    
    if nargin<4
        LowLim=[];
        HighLim=[];
    end
    
    var=struct('Indx',findcells(r,'FamName',familyname{nstru}),...
               'Parameter',{fieldtochange{nstru}},...
               'LowLim',LowLim,...
               'HighLim',HighLim...
               );
    
    variable=[variable, var]; %#ok<*AGROW>
end


return