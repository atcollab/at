function [newring,radelemIndex,cavitiesIndex,energy] = atenable_6d(ring,varargin)

% Process the input arguments
[allpass,varargs]=getoption(varargin,'allpass',[]);
[octupass,varargs]=getoption(varargs,'octupass',default_pass(''));
[sextupass,varargs]=getoption(varargs,'sextupass',default_pass(''));
[quadpass,varargs]=getoption(varargs,'quadpass',default_pass('auto'));
[wigglerpass,varargs]=getoption(varargs,'wigglerpass',default_pass('auto'));
[bendpass,varargs]=getoption(varargs,'bendpass',default_pass('auto'));
[cavipass,varargs]=getoption(varargs,'cavipass',default_pass('auto'));
[quantdiffpass,varargs]=getoption(varargs,'quantdiffpass',default_pass('auto'));
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
                elem.PassMethod=strrep(elem.PassMethod,'Pass','RadPass');
                elem.Energy=energy;
            end
            fprintf('%10s: %s\n', elem.FamName, elem.PassMethod);
        end
        
        function elem=fixelem(elem)
            % Explicit multipole modification
            elem.PassMethod=newpass;
            elem.Energy=energy;
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

    function defpass=default_pass(defpass)
        % Substitute the default pass method if 'allpass' is specified
        if ischar(allpass)
            defpass=allpass;
        end
    end
end

