function variable=atVariableBuilder(r,familyname,fieldtochange,LowLim,HighLim)
% 
% defines a simple variable structure for use with atmatch
% 
% set limits if not given and formats structure correctly.
% 
% familyname={'QD','SF','Dipole1',...}
% fieldtochange={{'PolynomB',{1,2}},{'PolynomB',{1,3}},{'BendingAngle',{1}},...}
% 
% if fieldtochange is a a cell array of one element, the same fieldtochange
% will be applied to all familynames.
% 
% familyname={@(r,v)function(r,v,...),'SF','Dipole1',...}
% fieldtochange={init_v_values,{'PolynomB',{1,3}},{'BendingAngle',{1}},...}
% 
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
% ex : v=atVariableBuilder(r,{'quad','sext','dip'},{{'PolynomB',{1,2}},{'PolynomB',{1,3}},{'BendingAngle'}})
% ex : v=atVariableBuilder(r,{'Q1','Q2','Q3'},{{'PolynomB',{1,2}}})
% ex : v=atVariableBuilder(r,familyname,fieldtochange,LowLim,HighLim)
% ex : v=atVariableBuilder(r,{@(r,~)matchchrom(r,fam1,fam2)},{[]})
%                                 

% history of changes
% created 21-03-2013
% update 29-03-2013 create many variables with the same parameter field.
% update 30-03-2013 create function variables.


variable=[];
if length(fieldtochange)==1
    fieldtochange=repmat(fieldtochange,1,length(familyname));
end
    
for nstru=1:length(familyname)
    if nargin<4
        LowLim=[];
        HighLim=[];
    end
    
    if isa(familyname{nstru},'function_handle') || isa(familyname{nstru},'numeric')
        var=struct('Indx',familyname{nstru},...
            'Parameter',{fieldtochange{nstru}},...
            'LowLim',LowLim,...
            'HighLim',HighLim...
            );
    else
        var=struct('Indx',findcells(r,'FamName',familyname{nstru}),...
            'Parameter',{fieldtochange{nstru}},...
            'LowLim',LowLim,...
            'HighLim',HighLim...
            );
        
    end

    variable=[variable, var]; %#ok<*AGROW>
end


return