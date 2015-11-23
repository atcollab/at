function variable=atVariableBuilder(varargin)
%atVarBuilder   create a simple variable structure for use with atmatch
%
% Single variable : it corresponds to a scalar numeric value to be varied in
% the optimization process. It may be applied to several elements.It is
% represented as a scalar structure.
%
%   var=atVariableBuilder(refpts,parameter,highlim,lowlim)
%       refpts:     indices of the variable elements or logical mask
%       parameter:	cell array defining the field name and indices of the
%                   variable parameter
%       lowlim:     minimum parameter value (default: no limit)
%       highlim:    maximum parameter value (default: no limit)
%
%       Example:	qf=atgetcells(ring,'FamName','QF');
%                   var=atVariableBuilder(qf,{'PolynomB',{2}});
%
%   var=atVariableBuilder(@func,inival,highlim,lowlim)
%       func:       function building a new ring for the given variable value
%                   called as new_ring=func(base_ring,variable)
%       inival:     initial value of the variable
%       lowlim:     minimum parameter value (default: no limit)
%       highlim:    maximum parameter value (default: no limit)
%
%       Example: var=atVariableBuilder(@(r,v) some_function(r,v,...), 0.0);
%
%   var=atVariableBuilder(ring,location,...)
%       In this syntax, the location may be specified as the family name of the
%       variable elements
%
%       Example: var=atVariableBuilder(ring,'QF',{'PolynomB',{2}});
%
% Multiple variables: if location,parameter,lowlim and highlim are cell arrays
% with the same length or with length 1, atVariableBuilder will build a
% structure array of variables. Examples:
%
%   vars=atVariableBuilder(ring,{'QD','SF'},{{'PolynomB',{1,2}},{'PolynomB',{1,3}}});
%
%   qf=atgetcells(ring,'FamName','QF');
%   qd=atgetcells(ring,'FamName','QD');
%   vars=atVariableBuilder({qf,qd},{{'PolynomB',{1,2}}});
%
%   vars=atVariableBuilder({qf,@buildring},{{'PolynomB',{1,2}},0.0})
%
% More sophisticated variables, can be defined using directly the variable
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

% history of changes
% created 21-03-2013
% update 29-03-2013 create many variables with the same parameter field.
% update 30-03-2013 create function variables.
% update 13-11-2015 reorganize function and help

if iscell(varargin{1}) && isfield(varargin{1}{1},'PassMethod')
    variable=atVariableBuilder(getid(varargin{1},varargin{2}),varargin{3:end});
elseif iscell(varargin{1})
    vals={{[]},{{}},{[]},{[]}};
    vals(1:nargin)=varargin;
    expand=1:max(cellfun(@length,varargin));
    location(expand)=cellfun(@tonum,vals{1},'UniformOutput',false);
    parameters(expand)=vals{2};
    lowl(expand)=vals{3};
    highl(expand)=vals{4};
    variable=struct('Indx',location,...
        'Parameter',parameters,...
        'LowLim',lowl,...
        'HighLim',highl...
        );
else
    vals=cellfun(@(i) {i}, varargin, 'UniformOutput',false);
    variable=atVariableBuilder(vals{:});
end

    function id=getid(ring,name)
        if iscell(name)
            id=cellfun(@(nm) getid(ring,nm), name, 'UniformOutput',false);
        elseif ischar(name)
            id=atgetcells(ring,'FamName',name);
        else
            id=name;
        end
    end
    function vnum=tonum(vnum)
        if islogical(vnum)
            vnum=find(vnum);
        end
    end

end
