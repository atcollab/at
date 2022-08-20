function [newring,radelemIndex,cavitiesIndex] = atdisable_6d(ring,varargin)

% Process the input arguments
[allpass,varargs]=getoption(varargin,'allpass',[]);
[octupass,varargs]=getoption(varargs,'octupass',default_pass('auto'));
[sextupass,varargs]=getoption(varargs,'sextupass',default_pass('auto'));
[quadpass,varargs]=getoption(varargs,'quadpass',default_pass('auto'));
[wigglerpass,varargs]=getoption(varargs,'wigglerpass',default_pass('auto'));
[bendpass,varargs]=getoption(varargs,'bendpass',default_pass('auto'));
[cavipass,varargs]=getoption(varargs,'cavipass',default_pass('auto'));
[quantdiffpass,varargs]=getoption(varargs,'quantdiffpass',default_pass('auto'));
[cavipass,bendpass,quadpass]=getargs(varargs,cavipass,bendpass,quadpass);

% Build the modification table
modfun.RFCavity=autoRFPass(cavipass);
modfun.Bend=autoMultipolePass(bendpass);
modfun.Quadrupole=autoMultipolePass(quadpass);
modfun.Sextupole=autoMultipolePass(sextupass);
modfun.Octupole=autoMultipolePass(octupass);
modfun.Wiggler=autoMultipolePass(wigglerpass);
modfun.QuantDiff=autoElemPass(quantdiffpass,'IdentityPass');
modfun.Other=@(elem) elem;

% Generate the new lattice
newring=cellfun(@modelem,ring,'UniformOutput',false);

if nargout > 1
    cavitiesIndex=atgetcells(newring,'PassMethod',@(elem,pass) endsWith(pass,'CavityPass'));
    radelemIndex=atgetcells(newring,'PassMethod',@(elem,pass) endsWith(pass,'RadPass'));
end

    function elem=modelem(elem)
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
            modfun=@fixelem;
        end
        
        function elem=varelem(elem)
            % 'auto' multipole modification
            elem.PassMethod=strrep(elem.PassMethod,'RadPass','Pass');
            if isfield(elem,'Energy'), elem=rmfield(elem,'Energy'); end
            fprintf('%10s: %s\n', elem.FamName, elem.PassMethod);
        end
        
        function elem=fixelem(elem)
            % Explicit multipole modification
            elem.PassMethod=newpass;
            if isfield(elem,'Energy'), elem=rmfield(elem,'Energy'); end
            fprintf('%10s: %s\n', elem.FamName, elem.PassMethod);
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
            fprintf('%10s: %s\n', elem.FamName, elem.PassMethod);
        end
    end

    function modfun=autoRFPass(newpass)
        % Returns the RF processing function
        if isempty(newpass)
            modfun=@(elem) elem;
        elseif strcmp(newpass, 'auto')
            modfun=@varelem;
        else
            modfun=@fixelem;
        end
        
        function elem=varelem(elem)
            % 'auto' RF modification
            if elem.Length > 0
                elem.PassMethod='DriftPass';
            else
                elem.PassMethod='IdentityPass';
            end
            fprintf('%10s: %s\n', elem.FamName, elem.PassMethod);
        end
        
        function elem=fixelem(elem)
            % Explicit RF modification
            elem.PassMethod=newpass;
            fprintf('%10s: %s\n', elem.FamName, elem.PassMethod);
        end
    end                

    function defpass=default_pass(defpass)
        % Substitute the default pass method if 'allpass' is specified
        if ischar(allpass)
            defpass=allpass;
        end
    end
end

