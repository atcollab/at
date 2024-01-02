function [newring,radelemIndex,cavitiesIndex] = atdisable_6d(ring,varargin)
%ATDISABLE_6D  switches radiation and cavity off
%
% [NEWRING,RADINDEX,CAVINDEX] = ATDISABLE_6D(RING,CAVIPASS,BENDPASS,QUADPASS)
%    Changes passmethods to turn off cavities, radiation damping and all
%    elements acting on the particle momentum.
%
% The default is to turn everything OFF., 
%
%  INPUTS:
%  1. RING      initial AT structure
%  2. CAVIPASS  pass method for cavities
%               '' makes no change,
%               'auto' sets'IdentityPass' or 'DriftPass' depending of cavity length)
%               anything else is used as the new PassMethod.
%  3. BENDPASS  pass method for bending magnets
%               '' makes no change,
%               'auto' substitutes 'RadPass' with 'Pass' in any method
%               anything else is used as the new PassMethod.
%  4. QUADPASS  pass method for quadrupoles
%               '' makes no change,
%               'auto' substitutes 'RadPass' with 'Pass' in any method
%               anything else is used as the new PassMethod.
%
%  [...] = ATDISABLE_6D(...[,keyword,value]...)
%   The following keywords trigger the processing of the following elements:
%
%   'allpass'        Defines the default pass method for all elements not
%                    explicitly specified. Replaces the following default
%                    values.
%   'cavipass'       pass method for RF cavities. Default 'auto'
%   'bendpass'       pass method for bending magnets. Default 'auto'
%   'quadpass'       pass method for quadrupoles. Default 'auto'
%   'sextupass'      pass method for sextupoles. Default 'auto'
%   'octupass'       pass method for bending magnets. Default 'auto'
%   'wigglerpass'	 pass method for wigglers. Default 'auto'
%   'quantdiffpass'  pass method for quantum diffusion. Default 'auto'
%   'energylosspass' pass method for atenergyloss element. Default 'auto'
%   'simplequantdiffpass' pass method for SimpleQuantDiff element. Default 'auto'
%   'simpleradiationpass' pass method for SimpleRadiation element. Default 'auto'
%
%   OUPUTS:
%   1. NEWRING   Output ring
%   2. RADINDEX  Indices of elements with radiation
%   3. CAVINDEX  Indices of active cavities
%
%  EXAMPLES:
%
%>> ringrad=atdisable_6d(ring);
%   Turns off all elements acting on momentum.
%
%>> ringrad=atdisable_6d(ring,'auto','allpass','');
%   Turns cavities off and leaves everything else unchanged.
%
%>> ringrad=atdisable_6d(ring,'allpass','auto','cavipass','');
%   Turns off everything except RF cavities.
%
%  See also ATENABLE_6D, CHECK_6D, ATCAVITYON, ATCAVITYOFF

% Process the keyword arguments
[allpass,varargs]=getoption(varargin,'allpass',[]);
[octupass,varargs]=getoption(varargs,'octupass',default_pass('auto'));
[sextupass,varargs]=getoption(varargs,'sextupass',default_pass('auto'));
[quadpass,varargs]=getoption(varargs,'quadpass',default_pass('auto'));
[wigglerpass,varargs]=getoption(varargs,'wigglerpass',default_pass('auto'));
[bendpass,varargs]=getoption(varargs,'bendpass',default_pass('auto'));
[cavipass,varargs]=getoption(varargs,'cavipass',default_pass('auto'));
[quantdiffpass,varargs]=getoption(varargs,'quantdiffpass',default_pass('auto'));
[energylosspass,varargs]=getoption(varargs,'energylosspass',default_pass('auto'));
[simplequantdiffpass,varargs]=getoption(varargs,'simplequantdiffpass',default_pass('auto'));
[simpleradiationpass,varargs]=getoption(varargs,'simpleradiationpass',default_pass('auto'));
% Process the positional arguments
[cavipass,bendpass,quadpass]=getargs(varargs,cavipass,bendpass,quadpass);

% Build the modification table
modfun.RFCavity=autoIdentityPass(cavipass);
modfun.Bend=autoMultipolePass(bendpass);
modfun.Quadrupole=autoMultipolePass(quadpass);
modfun.Sextupole=autoMultipolePass(sextupass);
modfun.Octupole=autoMultipolePass(octupass);
modfun.Wiggler=autoMultipolePass(wigglerpass);
modfun.QuantDiff=autoIdentityPass(quantdiffpass);
modfun.EnergyLoss=autoIdentityPass(energylosspass);
modfun.SimpleQuantDiff=autoIdentityPass(simplequantdiffpass);
modfun.SimpleRadiation=autoIdentityPass(simpleradiationpass);
modfun.Other=@(elem) elem;

% Generate the new lattice
newring=cellfun(@modelem,ring,'UniformOutput',false);

if nargout > 1
    cavitiesIndex=atgetcells(newring,'PassMethod',@(elem,pass) endsWith(pass,'CavityPass'));
    radelemIndex=atgetcells(newring,'PassMethod',@(elem,pass) endsWith(pass,'RadPass'));
end

    function elem=modelem(elem)
        %Modify the tracking PassMethod removing radiation
        cls=getclass_6d(elem);
        if any(strcmp(cls,fieldnames(modfun)))
            elem=modfun.(cls)(elem);
        end
    end

    function modfun=autoMultipolePass(newpass)
        % Returns the multipole processing function
        if isempty(newpass)
            modfun=@(elem) elem;
        elseif strcmp(newpass, 'auto')
            modfun=@varelem;
        else
            modfun=setpassenergy(newpass);
        end
        
        function elem=varelem(elem)
            % 'auto' multipole modification
            strrep(elem.PassMethod,'QuantPass','Pass)');
            elem.PassMethod=strrep(strrep(elem.PassMethod,'QuantPass','Pass'),'RadPass','Pass');
            if isfield(elem,'Energy'), elem=rmfield(elem,'Energy'); end
        end
    end

    function modfun=autoIdentityPass(newpass)
        % Returns the RF processing function
        if isempty(newpass)
            modfun=@(elem) elem;
        elseif strcmp(newpass, 'auto')
            modfun=@varelem;
        else
            modfun=setpass(newpass);
        end
        
        function elem=varelem(elem)
            % 'auto' RF modification
            if elem.Length > 0
                elem.PassMethod='DriftPass';
            else
                elem.PassMethod='IdentityPass';
            end
        end
    end

    function setfun=setpass(npass)
        function elem=newelem(elem)
            elem.PassMethod=npass;
        end
        setfun=@newelem;
    end

    function setfun=setpassenergy(npass)
        function elem=newelem(elem)
            elem.PassMethod=npass;
            if isfield(elem,'Energy'), elem=rmfield(elem,'Energy'); end
        end
        setfun=@newelem;
    end

    function defpass=default_pass(defpass)
        % Substitute the default pass method if 'allpass' is specified
        if ischar(allpass)
            defpass=allpass;
        end
    end
end

