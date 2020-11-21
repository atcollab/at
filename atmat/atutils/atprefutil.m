function varargout = atprefutil(action,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

persistent prefs

if ~isstruct(prefs)
    prefs=initoptions();
end

switch action
    case 'set'
        prset(varargin{:});
        varargout={};
    case 'get'
        varargout={prget(varargin{:})};
end

    function prset(pref,value)
        if nargin < 1               % Reset all
            prefs=initoptions();
        elseif isfield(prefs,pref)
            if nargin < 2           % Reset the selected constant
                iniv=initoptions();
                value=iniv.(pref);
            end
            prefs.(pref)=value;     % Set the selected constant
        else
            error('AT:SetPref','No %s constant is defined',pref);
        end
    end

    function value=prget(pref)
        if nargin < 1               % Get all
            value=prefs;
        else                        % Get the selected constant
            try
                value=prefs.(pref);
            catch
                error('AT:GetPref','No % constant is defined',pref);
            end
        end
    end
end
