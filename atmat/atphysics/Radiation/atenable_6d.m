function [newring,radelemIndex,cavitiesIndex,energy] = atenable_6d(ring,varargin)
%ATENABLE_6D switches RF and radiation on
%
%[NEWRING,RADINDEX,CAVINDEX] = ATENABLE_6D(RING,CAVIPASS,BENDPASS,QUADPASS)
%    Changes passmethods to get RF cavity acceleration and radiation
%    damping.
%
% The default is to turn cavities ON and set radiation in dipoles,
% quadrupoles and wigglers.
%
%  INPUTS:
%  1. RING      initial AT structure
%  2. CAVIPASS	pass method for cavities
%               '' makes no change,
%               'auto' set 'RFCavityPass',
%               anything else is used as the new PassMethod.
%  3. BENDPASS	pass method for bending magnets
%               '' makes no change,
%               'auto' substitutes 'Pass' with 'RadPass' in any method,
%               anything else is used as the new PassMethod.
%  4. QUADPASS	pass method for quadrupoles
%               '' makes no change,
%               'auto' substitutes 'Pass' with 'RadPass' in any method,
%               anything else is used as the new PassMethod.
%
%  [...] = ATENABLE_6D(...,keyword,value)
%   The following keywords trigger the processing of the following elements:
%
%   'allpass'       Defines the default pass method for all elements not
%                   explicitly specified. Replaces the following default
%                   values.
%   'cavipass'      pass method for RF cavities. Default 'auto'
%   'bendpass'      pass method for bending magnets. Default 'auto'
%   'quadpass'      pass method for quadrupoles. Default 'auto'
%   'sextupass'     pass method for sextupoles. Default ''
%   'octupass'      pass method for octupoles. Default ''
%   'wigglerpass'   pass method for wigglers. Default 'auto'
%   'quantdiffpass' pass method for quantum radiation. default 'auto'
%
%  OUPUTS:
%  1. NEWRING   Output ring
%  2. RADINDEX  Indices of elements with radiation
%  3. CAVINDEX  Indices of active cavities
%
%  EXAMPLES:
%
%>> ringrad=atenable_6d(ring);
%   Turns cavities on and sets radiation in bending magnets, quadrupoles and wigglers (default)
%
%>> ringrad=atenable_6d(ring,'auto','allpass','');
%   Turns cavities on and leaves everything else unchanged
%
%>> ringrad=atenable_6d(ring,'allpass','','bendpass','auto');
%   Turns on radiation in bending magnets and leaves everything else unchanged
%
%  See also ATDISABLE_6D, CHECK_6D, ATCAVITYON, ATCAVITYOFF

% Process the keyword arguments
[allpass,varargs]=getoption(varargin,'allpass',[]);
[octupass,varargs]=getoption(varargs,'octupass',default_pass(''));
[sextupass,varargs]=getoption(varargs,'sextupass',default_pass(''));
[quadpass,varargs]=getoption(varargs,'quadpass',default_pass('auto'));
[wigglerpass,varargs]=getoption(varargs,'wigglerpass',default_pass('auto'));
[bendpass,varargs]=getoption(varargs,'bendpass',default_pass('auto'));
[cavipass,varargs]=getoption(varargs,'cavipass',default_pass('auto'));
[quantdiffpass,varargs]=getoption(varargs,'quantdiffpass',default_pass('auto'));
% Process the positional arguments
[cavipass,bendpass,quadpass]=getargs(varargs,cavipass,bendpass,quadpass);

energy=atenergy(ring);

% Build the modification table
mod.RFCavity=autoElemPass(cavipass,'RFCavityPass');
mod.Bend=autoMultipolePass(bendpass);
mod.Quadrupole=autoMultipolePass(quadpass);
mod.Sextupole=autoMultipolePass(sextupass);
mod.Octupole=autoMultipolePass(octupass);
mod.Wiggler=autoMultipolePass(wigglerpass);
mod.QuantDiff=autoElemPass(quantdiffpass,'QuantDiffPass');
mod.Other=@(elem) elem;

% Generate the new lattice
newring=cellfun(@modelem,ring,'UniformOutput',false);

if nargout > 1
    cavitiesIndex=atgetcells(newring,'PassMethod',@(elem,pass) endsWith(pass,'CavityPass'));
    radelemIndex=atgetcells(newring,'PassMethod',@(elem,pass) endsWith(pass,'RadPass'));
end

    function elem=modelem(elem)
        cls=getclass_6d(elem);
        if any(strcmp(cls,fieldnames(mod)))
            elem=mod.(cls)(elem);
        end
    end

    function modfun=autoMultipolePass(newpass)
        % Returns the multipole processing function
        if isempty(newpass)
            modfun=@(elem) elem;
        elseif strcmp(newpass, 'auto')
            modfun=@varelem;
        else
            modfun=@fixelem;
        end
        
        function elem=varelem(elem)
            % 'auto' multipole modification
            if ~(endsWith(elem.PassMethod,'RadPass') || (isfield(elem,'BendingAngle') && elem.BendingAngle == 0))
                elem.PassMethod=strrep(strrep(elem.PassMethod,'QuantPass','Pass'),'Pass','RadPass');
                elem.Energy=energy;
            end
        end
        
        function elem=fixelem(elem)
            % Explicit multipole modification
            elem.PassMethod=newpass;
            elem.Energy=energy;
        end
    end

    function modfun=autoElemPass(newpass,defpass)
        % Returns the generic processing function
        if isempty(newpass)
            modfun=@(elem) elem;
        else
            if strcmp(newpass, 'auto')
                newpass=defpass;
            end
            modfun=@modelem;
        end
        
        function elem=modelem(elem)
            % Default element modification
            elem.PassMethod=newpass;
        end
    end

    function defpass=default_pass(defpass)
        % Substitute the default pass method if 'allpass' is specified
        if ischar(allpass)
            defpass=allpass;
        end
    end
end

